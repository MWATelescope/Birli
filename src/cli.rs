//! Command Line Interface helpers for Birli

use crate::{
    error::{BirliError, BirliError::DryRun, CLIError::InvalidCommandLineArgument},
    flags::FlagContext,
    io::{aocal::AOCalSols, IOContext},
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        hifitime::Epoch,
        io::{ms::MeasurementSetWriter, uvfits::UvfitsWriter, VisWritable},
        mwalib::{CorrelatorContext, GeometricDelaysApplied},
        ndarray::s,
        precession::{precess_time, PrecessionInfo},
        LatLngHeight, RADec,
    },
    passband_gains::{PFB_COTTER_2014_10KHZ, PFB_JAKE_2022_200HZ},
    with_increment_duration, Axis, Complex, FlagFileSet, PreprocessContext, VisSelection,
};
use cfg_if::cfg_if;
use clap::{arg, command, ErrorKind::ArgumentNotFound, PossibleValue, ValueHint::FilePath};
use itertools::{izip, Itertools};
use log::{debug, info, trace, warn};
use marlu::{Jones, MwaObsContext, ObsContext, VisContext};
use prettytable::{cell, format as prettyformat, row, table};
use std::{
    collections::HashMap,
    convert::Into,
    env,
    ffi::OsString,
    fmt::{Debug, Display},
    time::Duration,
};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use aoflagger_sys::{cxx_aoflagger_new};
    }
}

/// Args for preprocessing a correlator context.
pub struct BirliContext {
    /// mwalib::CorrelatorContext
    pub corr_ctx: CorrelatorContext,
    /// Preprocessing parameters
    pub prep_ctx: PreprocessContext,
    /// selected visibility indices
    pub vis_sel: VisSelection,
    /// Flagging Parameters
    pub flag_ctx: FlagContext,
    /// Input / output paths
    pub io_ctx: IOContext,
    /// temporal averaging factor
    pub avg_time: usize,
    /// spectral averaging factor
    pub avg_freq: usize,
    /// temporal chunking factor
    pub num_timesteps_per_chunk: Option<usize>,
}

// Add build-time information from the "built" crate.
include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// stolen from hyperdrive
/// Write many info-level log lines of how this executable was compiled.
///
/// # Errors
///
/// propagates writeln! fails
pub fn fmt_build_info(f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match GIT_HEAD_REF {
        Some(hr) => {
            let dirty = GIT_DIRTY.unwrap_or(false);
            writeln!(
                f,
                "Compiled on git commit hash: {}{}",
                GIT_COMMIT_HASH.unwrap(),
                if dirty { " (dirty)" } else { "" }
            )?;
            writeln!(f, "            git head ref: {}", hr)?;
        }
        None => writeln!(f, "Compiled on git commit hash: <no git info>")?,
    }
    writeln!(f, "            {}", BUILT_TIME_UTC)?;
    writeln!(f, "         with compiler {}", RUSTC_VERSION)?;
    writeln!(f)?;
    Ok(())
}

fn time_details(
    gps_time_ms: u64,
    phase_centre: RADec,
    array_pos: LatLngHeight,
) -> (String, String, f64, PrecessionInfo) {
    let epoch = Epoch::from_gpst_seconds(gps_time_ms as f64 / 1e3);
    let (y, mo, d, h, mi, s, ms) = epoch.as_gregorian_utc();
    let precession_info = precess_time(
        phase_centre,
        epoch,
        array_pos.longitude_rad,
        array_pos.latitude_rad,
    );
    (
        format!("{:02}-{:02}-{:02}", y, mo, d),
        format!(
            "{:02}:{:02}:{:02}.{:03}",
            h,
            mi,
            s,
            (ms as f64 / 1e6).round()
        ),
        epoch.as_mjd_utc_seconds(),
        precession_info,
    )
}

impl Display for BirliContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{} version {}",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION"),
        )?;

        fmt_build_info(f)?;

        writeln!(
            f,
            "observation name:     {}",
            self.corr_ctx.metafits_context.obs_name
        )?;

        writeln!(f, "Array position:       {}", &self.prep_ctx.array_pos)?;
        writeln!(f, "Phase centre:         {}", &self.prep_ctx.phase_centre)?;
        let pointing_centre = RADec::from_mwalib_tile_pointing(&self.corr_ctx.metafits_context);
        if pointing_centre != self.prep_ctx.phase_centre {
            writeln!(f, "Pointing centre:      {}", &pointing_centre)?;
        }

        let coarse_chan_flag_idxs: Vec<usize> = self
            .flag_ctx
            .coarse_chan_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        // TODO: actually display this.
        let _fine_chan_flag_idxs: Vec<usize> = self
            .flag_ctx
            .fine_chan_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        let timestep_flag_idxs: Vec<usize> = self
            .flag_ctx
            .timestep_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();
        let ant_pairs = self.vis_sel.get_ant_pairs(&self.corr_ctx.metafits_context);
        #[allow(clippy::needless_collect)]
        let baseline_flag_idxs: Vec<usize> = self
            .flag_ctx
            .get_baseline_flags(&ant_pairs)
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .collect();

        let (sched_start_date, sched_start_time, sched_start_mjd_s, sched_start_prec) =
            time_details(
                self.corr_ctx.metafits_context.sched_start_gps_time_ms,
                self.prep_ctx.phase_centre,
                self.prep_ctx.array_pos,
            );
        writeln!(
            f,
            "Scheduled start:      {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
            sched_start_date, sched_start_time,
            self.corr_ctx.metafits_context.sched_start_unix_time_ms as f64 / 1e3,
            self.corr_ctx.metafits_context.sched_start_gps_time_ms as f64 / 1e3,
            sched_start_mjd_s,
            sched_start_prec.lmst.to_degrees(),
            sched_start_prec.lmst_j2000.to_degrees(),
            sched_start_prec.array_latitude_j2000.to_degrees(),
        )?;
        let (sched_end_date, sched_end_time, sched_end_mjd_s, sched_end_prec) = time_details(
            self.corr_ctx.metafits_context.sched_end_gps_time_ms,
            self.prep_ctx.phase_centre,
            self.prep_ctx.array_pos,
        );
        writeln!(
            f,

            "Scheduled end:        {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
            sched_end_date, sched_end_time,
            self.corr_ctx.metafits_context.sched_end_unix_time_ms as f64 / 1e3,
            self.corr_ctx.metafits_context.sched_end_gps_time_ms as f64 / 1e3,
            sched_end_mjd_s,
            sched_end_prec.lmst.to_degrees(),
            sched_end_prec.lmst_j2000.to_degrees(),
            sched_end_prec.array_latitude_j2000.to_degrees(),
        )?;
        let int_time_s = self.corr_ctx.metafits_context.corr_int_time_ms as f64 / 1e3;
        let sched_duration_s = self.corr_ctx.metafits_context.sched_duration_ms as f64 / 1e3;
        writeln!(
            f,
            "Scheduled duration:   {:.3}s = {:3} * {:.3}s",
            sched_duration_s,
            (sched_duration_s / int_time_s).ceil(),
            int_time_s
        )?;
        let quack_duration_s = self.corr_ctx.metafits_context.quack_time_duration_ms as f64 / 1e3;
        writeln!(
            f,
            "Quack duration:       {:.3}s = {:3} * {:.3}s",
            quack_duration_s,
            (quack_duration_s / int_time_s).ceil(),
            int_time_s
        )?;
        let num_avg_timesteps =
            (self.vis_sel.timestep_range.len() as f64 / self.avg_time as f64).ceil() as usize;
        let avg_int_time_s = int_time_s * self.avg_time as f64;
        writeln!(
            f,
            "Output duration:      {:.3}s = {:3} * {:.3}s{}",
            num_avg_timesteps as f64 * avg_int_time_s,
            num_avg_timesteps,
            avg_int_time_s,
            if self.avg_time == 1 {
                "".into()
            } else {
                format!(" ({}x)", self.avg_time)
            }
        )?;

        let total_bandwidth_mhz = self.corr_ctx.metafits_context.obs_bandwidth_hz as f64 / 1e6;
        let fine_chan_width_khz =
            self.corr_ctx.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;
        let fine_chans_per_coarse = self
            .corr_ctx
            .metafits_context
            .num_corr_fine_chans_per_coarse;

        writeln!(
            f,
            "Scheduled Bandwidth:  {:.3}MHz = {:3} * {:3} * {:.3}kHz",
            total_bandwidth_mhz,
            self.corr_ctx.metafits_context.num_metafits_coarse_chans,
            fine_chans_per_coarse,
            fine_chan_width_khz
        )?;

        let out_bandwidth_mhz = self.vis_sel.coarse_chan_range.len() as f64
            * fine_chans_per_coarse as f64
            * fine_chan_width_khz
            / 1e3;
        let num_avg_chans = (self.vis_sel.coarse_chan_range.len() as f64
            * fine_chans_per_coarse as f64
            / self.avg_freq as f64)
            .ceil() as usize;
        let avg_fine_chan_width_khz = fine_chan_width_khz * self.avg_freq as f64;
        writeln!(
            f,
            "Output Bandwidth:     {:.3}MHz = {:9} * {:.3}kHz{}",
            out_bandwidth_mhz,
            num_avg_chans,
            avg_fine_chan_width_khz,
            if self.avg_freq == 1 {
                "".into()
            } else {
                format!(" ({}x)", self.avg_freq)
            }
        )?;

        let first_epoch =
            Epoch::from_gpst_seconds(self.corr_ctx.timesteps[0].gps_time_ms as f64 / 1e3);
        let (y, mo, d, ..) = first_epoch.as_gregorian_utc();

        let mut timestep_table = table!([
            "",
            format!("{:02}-{:02}-{:02} UTC +", y, mo, d),
            "unix [s]",
            "gps [s]",
            "p",
            "c",
            "g",
            "s",
            "f"
        ]);
        timestep_table.set_format(*prettyformat::consts::FORMAT_CLEAN);

        let provided_timestep_indices = self.corr_ctx.provided_timestep_indices.clone();
        let common_timestep_indices = self.corr_ctx.common_timestep_indices.clone();
        let common_good_timestep_indices = self.corr_ctx.common_good_timestep_indices.clone();
        for (timestep_idx, timestep) in self.corr_ctx.timesteps.iter().enumerate() {
            let provided = provided_timestep_indices.contains(&timestep_idx);
            let selected = self.vis_sel.timestep_range.contains(&timestep_idx);
            let common = common_timestep_indices.contains(&timestep_idx);
            let good = common_good_timestep_indices.contains(&timestep_idx);
            let flagged = timestep_flag_idxs.contains(&timestep_idx);

            let (_, time, ..) = time_details(
                timestep.gps_time_ms,
                self.prep_ctx.phase_centre,
                self.prep_ctx.array_pos,
            );
            let row = row![r =>
                format!("ts{}:", timestep_idx),
                time,
                format!("{:.3}", timestep.unix_time_ms as f64 / 1e3),
                format!("{:.3}", timestep.gps_time_ms as f64 / 1e3),
                if provided {"p"} else {""},
                if common {"c"} else {""},
                if good {"g"} else {""},
                if selected {"s"} else {""},
                if flagged {"f"} else {""}
            ];
            timestep_table.add_row(row);
        }

        writeln!(
            f,
            "Timestep details (all={}, provided={}, common={}, good={}, select={}, flag={}):\n{}",
            self.corr_ctx.num_timesteps,
            self.corr_ctx.num_provided_timesteps,
            self.corr_ctx.num_common_timesteps,
            self.corr_ctx.num_common_good_timesteps,
            self.vis_sel.timestep_range.len(),
            timestep_flag_idxs.len(),
            timestep_table
        )?;

        let mut coarse_chan_table = table!([
            "",
            "gpu",
            "corr",
            "rec",
            "cen [MHz]",
            "p",
            "c",
            "g",
            "s",
            "f"
        ]);
        coarse_chan_table.set_format(*prettyformat::consts::FORMAT_CLEAN);
        // coarse_chan_table
        let provided_coarse_chan_indices = self.corr_ctx.provided_coarse_chan_indices.clone();
        let common_coarse_chan_indices = self.corr_ctx.common_coarse_chan_indices.clone();
        let common_good_coarse_chan_indices = self.corr_ctx.common_good_coarse_chan_indices.clone();
        for (chan_idx, chan) in self.corr_ctx.coarse_chans.iter().enumerate() {
            let provided = provided_coarse_chan_indices.contains(&chan_idx);
            let selected = self.vis_sel.coarse_chan_range.contains(&chan_idx);
            let common = common_coarse_chan_indices.contains(&chan_idx);
            let good = common_good_coarse_chan_indices.contains(&chan_idx);
            let flagged = coarse_chan_flag_idxs.contains(&chan_idx);
            let row = row![r =>
                format!("cc{}:", chan_idx),
                chan.gpubox_number,
                chan.corr_chan_number,
                chan.rec_chan_number,
                format!("{:.4}", chan.chan_centre_hz as f64 / 1e6),
                if provided {"p"} else {""},
                if common {"c"} else {""},
                if good {"g"} else {""},
                if selected {"s"} else {""},
                if flagged {"f"} else {""}
            ];
            coarse_chan_table.add_row(row);
        }

        writeln!(
            f,
            "Coarse channel details (metafits={}, provided={}, common={}, good={}, select={}, flag={}):\n{}",
            self.corr_ctx.num_coarse_chans,
            self.corr_ctx.num_provided_coarse_chans,
            self.corr_ctx.num_common_coarse_chans,
            self.corr_ctx.num_common_good_coarse_chans,
            self.vis_sel.coarse_chan_range.len(),
            coarse_chan_flag_idxs.len(),
            coarse_chan_table
        )?;

        writeln!(
            f,
            "Antenna details (all={}, flag={})",
            self.corr_ctx.metafits_context.num_ants,
            self.flag_ctx
                .antenna_flags
                .iter()
                .enumerate()
                .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
                .count(),
            // format!("\n{}", ant_table)
        )?;

        writeln!(
            f,
            "Baseline Details (all={}, auto={}, select={}, flag={}):",
            self.corr_ctx.metafits_context.num_baselines,
            self.corr_ctx.metafits_context.num_ants,
            self.vis_sel.baseline_idxs.len(),
            baseline_flag_idxs.len(),
        )?;

        // TODO: show free memory with https://docs.rs/sys-info/latest/sys_info/fn.mem_info.html

        let num_sel_timesteps = self.vis_sel.timestep_range.len();
        let num_sel_chans = self.vis_sel.coarse_chan_range.len() * fine_chans_per_coarse;
        let num_sel_baselines = self.vis_sel.baseline_idxs.len();
        let num_sel_pols = self.corr_ctx.metafits_context.num_visibility_pols;
        let mem_selected_bytes = self.vis_sel.estimate_bytes_best(fine_chans_per_coarse);
        let mem_per_timestep_gib =
            mem_selected_bytes as f64 / num_sel_timesteps as f64 / 1024.0_f64.powi(3);

        writeln!(
            f,
            "Estimated memory usage per timestep =           {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_chans,
            num_sel_baselines,
            std::mem::size_of::<Jones<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib,
        )?;

        if let Some(num_timesteps) = self.num_timesteps_per_chunk {
            writeln!(
                f,
                "Estimated memory per chunk          = {:5}ts * {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
                num_timesteps,
                num_sel_chans,
                num_sel_baselines,
                std::mem::size_of::<Jones<f32>>(),
                std::mem::size_of::<f32>(),
                std::mem::size_of::<bool>(),
                mem_per_timestep_gib * num_timesteps as f64,
            )?;
        }

        writeln!(
            f,
            "Estimated memory selected           = {:5}ts * {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_timesteps,
            num_sel_chans,
            num_sel_baselines,
            std::mem::size_of::<Jones<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib * num_sel_timesteps as f64,
        )?;

        let avg_mem_per_timestep_gib = (num_avg_chans
            * num_sel_baselines
            * num_sel_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())) as f64
            / 1024.0_f64.powi(3);

        writeln!(
            f,
            "Estimated output size               = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_avg_timesteps,
            num_avg_chans,
            num_sel_baselines,
            num_sel_pols,
            std::mem::size_of::<Complex<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            avg_mem_per_timestep_gib * num_avg_timesteps as f64,
        )?;

        writeln!(f, "Preprocessing Context: \n{}", &self.prep_ctx)?;

        Ok(())
    }
}

impl BirliContext {
    // TODO: try struct instead of builder
    #[allow(clippy::cognitive_complexity)]
    fn get_matches<I, T>(args: I) -> Result<clap::ArgMatches, BirliError>
    where
        I: IntoIterator<Item = T> + Debug,
        T: Into<OsString> + Clone,
    {
        #[allow(unused_mut)]
        let mut app = command!()
            .subcommand_precedence_over_arg(true)
            .arg_required_else_help(true)
            .next_line_help(false)
            .about("Preprocess Murchison Widefield Array MetaFITS and GPUFITS data \
                    into usable astronomy formats.")
            .args(&[
                // input options
                arg!(-m --metafits <PATH> "Metadata file for the observation")
                    .required(true)
                    .value_hint(FilePath)
                    .help_heading("INPUT"),
                arg!(fits_paths: <PATHS>... "GPUBox files to process")
                    .help_heading("INPUT")
                    .value_hint(FilePath)
                    .required(true),

                // processing options
                arg!(--"phase-centre" "Override Phase centre from metafits (degrees)")
                    .value_names(&["RA", "DEC"])
                    .required(false),
                arg!(--"pointing-centre" "Use pointing instead phase centre")
                    .conflicts_with("phase-centre"),
                arg!(--"emulate-cotter" "Use Cotter's array position, not MWAlib's"),
                arg!(--"dry-run" "Just print the summary and exit"),
                arg!(--"no-draw-progress" "do not show progress bars"),

                // selection options
                // TODO: make this work the same way as rust ranges. start <= x < end
                arg!(--"sel-time" "Timestep index range (inclusive) to select")
                    .help_heading("SELECTION")
                    .value_names(&["MIN", "MAX"])
                    .required(false),
                arg!(--"sel-ants" <ANTS>... "[WIP] Antenna to select")
                    .help_heading("SELECTION")
                    .multiple_values(true)
                    .required(false),
                arg!(--"no-sel-flagged-ants" "[WIP] Deselect flagged antennas")
                    .help_heading("SELECTION"),
                arg!(--"no-sel-autos" "[WIP] Deselect autocorrelations")
                    .help_heading("SELECTION"),

                // resource limit options
                arg!(--"time-chunk" <STEPS> "[WIP] Process observation in chunks of <STEPS> timesteps.")
                    .help_heading("RESOURCE LIMITS")
                    .required(false)
                    .conflicts_with("max-memory"),
                arg!(--"max-memory" <GIBIBYTES> "[WIP] Estimate --time-chunk with <GIBIBYTES> GiB each chunk.")
                    .help_heading("RESOURCE LIMITS")
                    .required(false),

                // flagging options
                // -> timesteps
                arg!(--"flag-init" <SECONDS> "[WIP] Flag <SECONDS> after first common time (quack time)")
                    .alias("--quack-time")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-init-steps" <COUNT> "[WIP] Flag <COUNT> steps after first common time")
                    .help_heading("FLAGGING")
                    .required(false)
                    .conflicts_with("flag-init"),
                arg!(--"flag-end" <SECONDS> "[WIP] Flag seconds before the last provided time")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-end-steps" <COUNT> "[WIP] Flag <COUNT> steps before the last provided")
                    .help_heading("FLAGGING")
                    .required(false)
                    .conflicts_with("flag-end"),
                arg!(--"flag-times" <STEPS>... "[WIP] Flag additional time steps")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                // -> channels
                arg!(--"flag-coarse-chans" <CHANS> ... "[WIP] Flag additional coarse chan indices")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                arg!(--"flag-edge-width" <KHZ> "[WIP] Flag bandwidth [kHz] at the ends of each coarse chan")
                    .help_heading("FLAGGING")
                    .required(false),
                arg!(--"flag-edge-chans" <COUNT> "[WIP] Flag <COUNT> fine chans on the ends of each coarse")
                    .help_heading("FLAGGING")
                    .conflicts_with("flag-edge-width")
                    .required(false),
                arg!(--"flag-fine-chans" <CHANS>... "[WIP] Flag fine chan indices in each coarse chan")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                arg!(--"flag-dc" "[WIP] Force flagging of DC centre chans")
                    .help_heading("FLAGGING"),
                arg!(--"no-flag-dc" "[WIP] Do not flag DC centre chans")
                    .help_heading("FLAGGING"),
                // -> antennas
                arg!(--"no-flag-metafits" "[WIP] Ignore antenna flags in metafits")
                    .help_heading("FLAGGING"),
                arg!(--"flag-antennas" <ANTS>... "[WIP] Flag antenna indices")
                    .help_heading("FLAGGING")
                    .multiple_values(true)
                    .required(false),
                // -> baselines
                arg!(--"flag-autos" "[WIP] Flag auto correlations")
                    .help_heading("FLAGGING"),

                // corrections
                arg!(--"no-cable-delay" "Do not perform cable length corrections")
                    .help_heading("CORRECTION"),
                arg!(--"no-geometric-delay" "Do not perform geometric corrections")
                    .help_heading("CORRECTION")
                    .alias("no-geom"),
                arg!(--"no-digital-gains" "Do not perform digital gains corrections")
                    .help_heading("CORRECTION"),
                arg!(--"passband-gains" <TYPE> "Type of PFB passband filter gains correction to apply")
                    .required(false)
                    .possible_values([
                        PossibleValue::new("none").help("No passband gains correction (unitary)"),
                        PossibleValue::new("cotter")
                            .help(
                                "_sb128ChannelSubbandValue2014FromMemo from
                                    subbandpassband.cpp in Cotter. Can only be used with resolutions of
                                    n * 10kHz"
                            ),
                        PossibleValue::new("jake")
                            .help("see: PFB_JAKE_2022_200HZ in src/passband_gains.rs"),
                    ])
                    .default_value("jake")
                    .alias("pfb-gains")
                    .help_heading("CORRECTION"),

                // calibration
                arg!(--"apply-di-cal" <PATH> "Apply DI calibration solutions before averaging")
                    .required(false)
                    .value_hint(FilePath),

                // averaging
                arg!(--"avg-time-res" <SECONDS> "Time resolution of averaged data")
                    .help_heading("AVERAGING")
                    .required(false),
                arg!(--"avg-time-factor" <FACTOR> "Average <FACTOR> timesteps per averaged timestep")
                    .help_heading("AVERAGING")
                    .required(false)
                    .conflicts_with("avg-time-res"),
                arg!(--"avg-freq-res" <KHZ> "Frequency resolution of averaged data")
                    .help_heading("AVERAGING")
                    .required(false),
                arg!(--"avg-freq-factor" <FACTOR> "Average <FACTOR> channels per averaged channel")
                    .help_heading("AVERAGING")
                    .required(false)
                    .conflicts_with("avg-freq-res"),

                // output options
                arg!(-f --"flag-template" <TEMPLATE> "The template used to name flag files. \
                        Percents are substituted for the zero-prefixed GPUBox ID, which can be up to \
                        3 characters long. Example: FlagFile%%%.mwaf")
                    .help_heading("OUTPUT")
                    .required(false),
                arg!(-u --"uvfits-out" <PATH> "Path for uvfits output")
                    .help_heading("OUTPUT")
                    .required(false),
                arg!(-M --"ms-out" <PATH> "Path for measurement set output")
                    .help_heading("OUTPUT")
                    .required(false),
            ]);
        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                app = app.args(&[
                    arg!(--"no-rfi" "Do not perform RFI Flagging with aoflagger")
                        .help_heading("AOFLAGGER"),
                    arg!(--"aoflagger-strategy" <PATH> "Strategy to use for RFI Flagging")
                        .value_hint(FilePath)
                        .help_heading("AOFLAGGER")
                        .required(false)
                ]);
            }
        };
        let matches = app.try_get_matches_from_mut(args)?;
        Ok(matches)
    }

    fn parse_io_matches(matches: &clap::ArgMatches) -> IOContext {
        IOContext {
            metafits_in: match matches.value_of_t("metafits") {
                Ok(path) => path,
                _ => unreachable!("--metafits <PATH> is required, enforced by clap"),
            },
            gpufits_in: match matches.values_of_t("fits_paths") {
                Ok(path) => path,
                _ => unreachable!("<PATHS> is required, enforced by clap"),
            },
            aocalsols_in: matches.value_of("apply-di-cal").map(Into::into),
            uvfits_out: matches.value_of("uvfits-out").map(Into::into),
            ms_out: matches.value_of("ms-out").map(Into::into),
            flag_template: matches.value_of("flag-template").map(Into::into),
        }
    }

    fn parse_vis_sel_matches(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
    ) -> Result<VisSelection, BirliError> {
        let mut vis_sel = VisSelection::from_mwalib(corr_ctx).unwrap();
        match matches
            .values_of_t::<usize>("sel-time")
            .map(|v| (v[0], v[1]))
        {
            Ok((from, to)) => {
                if from > to || to >= corr_ctx.num_timesteps {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--sel-time <FROM> <TO>".into(),
                        expected: format!("from <= to < num_timesteps={}", corr_ctx.num_timesteps),
                        received: format!("from={} to={}", from, to),
                    }));
                }
                vis_sel.timestep_range = from..(to + 1);
            }
            Err(err) => match err.kind() {
                ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        }
        Ok(vis_sel)
    }

    fn parse_flag_matches(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
    ) -> Result<FlagContext, BirliError> {
        let mut flag_ctx = FlagContext::from_mwalib(corr_ctx);
        match matches.values_of_t::<usize>("flag-times") {
            Ok(timestep_idxs) => {
                for (value_idx, &timestep_idx) in timestep_idxs.iter().enumerate() {
                    if timestep_idx >= corr_ctx.num_timesteps {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-times <TIMESTEPS>...".into(),
                            expected: format!(
                                "timestep_idx < num_timesteps={}",
                                corr_ctx.num_timesteps
                            ),
                            received: format!(
                                "timestep_idxs[{}]={}. all:{:?}",
                                value_idx, timestep_idx, timestep_idxs
                            ),
                        }));
                    }
                    flag_ctx.timestep_flags[timestep_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };
        match matches.values_of_t::<usize>("flag-coarse-chans") {
            Ok(coarse_chan_idxs) => {
                for (value_idx, &coarse_chan_idx) in coarse_chan_idxs.iter().enumerate() {
                    if coarse_chan_idx >= corr_ctx.num_coarse_chans {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-coarse-chans <CHANS>...".into(),
                            expected: format!(
                                "coarse_chan_idx < num_coarse_chans={}",
                                corr_ctx.num_coarse_chans
                            ),
                            received: format!(
                                "coarse_chan_idxs[{}]={}. all:{:?}",
                                value_idx, coarse_chan_idx, coarse_chan_idxs
                            ),
                        }));
                    }
                    flag_ctx.coarse_chan_flags[coarse_chan_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        match matches.values_of_t::<usize>("flag-fine-chans") {
            Ok(fine_chan_idxs) => {
                for (value_idx, &fine_chan_idx) in fine_chan_idxs.iter().enumerate() {
                    if fine_chan_idx >= fine_chans_per_coarse {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-fine-chans <CHANS>...".into(),
                            expected: format!(
                                "fine_chan_idx < num_fine_chans={}",
                                fine_chans_per_coarse
                            ),
                            received: format!(
                                "fine_chan_idxs[{}]={}. all:{:?}",
                                value_idx, fine_chan_idx, fine_chan_idxs
                            ),
                        }));
                    }
                    flag_ctx.fine_chan_flags[fine_chan_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };
        if matches.is_present("no-flag-metafits") {
            info!("Ignoring antenna flags from metafits.");
            // set antenna flags to all false
            flag_ctx.antenna_flags = vec![false; flag_ctx.antenna_flags.len()];
        }
        match matches.values_of_t::<usize>("flag-antennas") {
            Ok(antenna_idxs) => {
                for (value_idx, &antenna_idx) in antenna_idxs.iter().enumerate() {
                    if antenna_idx >= corr_ctx.metafits_context.num_ants {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-antennas <ANTS>...".into(),
                            expected: format!(
                                "antenna_idx < num_ants={}",
                                corr_ctx.metafits_context.num_ants
                            ),
                            received: format!(
                                "antenna_idxs[{}]={}. all:{:?}",
                                value_idx, antenna_idx, antenna_idxs
                            ),
                        }));
                    }
                    flag_ctx.antenna_flags[antenna_idx] = true;
                }
            }
            Err(err) => match err.kind() {
                ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };
        if matches.is_present("flag-autos") {
            flag_ctx.autos = true;
        }
        Ok(flag_ctx)
    }

    fn parse_avg_matches(
        matches: &clap::ArgMatches,
        corr_ctx: &CorrelatorContext,
    ) -> Result<(usize, usize), BirliError> {
        let avg_time: usize = match (
            matches.value_of_t::<usize>("avg-time-factor"),
            matches.value_of_t::<f64>("avg-time-res"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (_, Err(err)) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (Ok(_), Ok(_)) => {
                unreachable!("--avg-time-res conflicts with --avg-time-factor, enforced by clap")
            }
            (Ok(factor), _) => {
                if factor == 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-time-factor <FACTOR>".into(),
                        expected: "a positive, non-zero integer".into(),
                        received: format!("{}", factor),
                    }));
                }
                factor
            }
            (_, Ok(res)) => {
                let int_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1e3;
                let ratio = res / int_time_s;
                if ratio.is_infinite() || ratio.fract() > 1e-6 || ratio < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-time-res <RES>".into(),
                        expected: format!("a multiple of the integration time, {} [s]", int_time_s),
                        received: format!("{}", res),
                    }));
                }
                ratio.round() as _
            }
            _ => 1,
        };
        let avg_freq: usize = match (
            matches.value_of_t::<usize>("avg-freq-factor"),
            matches.value_of_t::<f64>("avg-freq-res"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (_, Err(err)) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (Ok(_), Ok(_)) => {
                unreachable!("--avg-freq-res conflicts with --avg-freq-factor, enforced by clap")
            }
            (Ok(factor), _) => {
                if factor == 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-freq-factor <FACTOR>".into(),
                        expected: "a positive, non-zero integer".into(),
                        received: format!("{}", factor),
                    }));
                }
                factor
            }
            (_, Ok(res)) => {
                let fine_chan_width_khz =
                    corr_ctx.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;
                let ratio = res / fine_chan_width_khz;
                if ratio.is_infinite() || ratio.fract() > 1e-6 || ratio < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--avg-freq-res <RES>".into(),
                        expected: format!(
                            "a multiple of the fine channel width, {} [KHz]",
                            fine_chan_width_khz
                        ),
                        received: format!("{}", res),
                    }));
                }
                ratio.round() as _
            }
            _ => 1,
        };
        Ok((avg_time, avg_freq))
    }

    fn parse_chunk_matches(
        corr_ctx: &CorrelatorContext,
        matches: &clap::ArgMatches,
        avg_time: usize,
        vis_sel: &VisSelection,
    ) -> Result<Option<usize>, BirliError> {
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let num_timesteps_per_chunk: Option<usize> = match (
            matches.value_of_t::<usize>("time-chunk"),
            matches.value_of_t::<f64>("max-memory"),
        ) {
            // filter any errors other than ArgumentNotFound
            (Err(err), _) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (_, Err(err)) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (Ok(_), Ok(_)) => {
                unreachable!("--time-chunk conflicts with --max-memory, enforced by clap")
            }
            (Ok(steps), _) => {
                if steps % avg_time != 0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--time-chunk <STEPS>".into(),
                        expected: format!(
                            "a multiple of the temporal averaging factor, {}",
                            avg_time
                        ),
                        received: format!("{}", steps),
                    }));
                }
                Some(steps)
            }
            (_, Ok(mem_gib)) => {
                let max_mem_bytes = mem_gib * 1024.0_f64.powi(3);
                if max_mem_bytes < 1.0 {
                    return Err(BirliError::CLIError(InvalidCommandLineArgument {
                        option: "--max-memory <GIBIBYTES>".into(),
                        expected: "at least one Byte".into(),
                        received: format!("{}B", max_mem_bytes),
                    }));
                }
                let bytes_selected = vis_sel.estimate_bytes_best(fine_chans_per_coarse);
                let bytes_per_timestep = bytes_selected / vis_sel.timestep_range.len();
                let bytes_per_avg_time = bytes_per_timestep * avg_time;
                if max_mem_bytes < bytes_selected as f64 {
                    if max_mem_bytes < bytes_per_avg_time as f64 {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--max-memory <GIBIBYTES>".into(),
                            expected: format!("at least enough memory for an averaged timestep ({} * {:.02} = {:.02} GiB)", avg_time, bytes_per_timestep as f64 / 1024.0_f64.powi(3), bytes_per_avg_time as f64 / 1024.0_f64.powi(3)),
                            received: format!("{}GiB", max_mem_bytes as f64 / 1024.0_f64.powi(3)),
                        }));
                    }
                    Some((max_mem_bytes / bytes_per_avg_time as f64).floor() as usize * avg_time)
                } else {
                    None
                }
            }
            _ => None,
        };

        // validate chunk size
        if let Some(chunk_size) = num_timesteps_per_chunk {
            assert!(
                matches.value_of("flag-template").is_none(),
                "chunking is not supported when writing .mwaf files using --flag-template"
            );
            info!("chunking output to {} timesteps per chunk", chunk_size);
        }

        Ok(num_timesteps_per_chunk)
    }

    fn parse_prep_matches(
        matches: &clap::ArgMatches,
        corr_ctx: &CorrelatorContext,
    ) -> Result<PreprocessContext, BirliError> {
        let mut prep_ctx = PreprocessContext {
            draw_progress: !matches.is_present("no-draw-progress"),
            ..PreprocessContext::default()
        };
        prep_ctx.array_pos = if matches.is_present("emulate-cotter") {
            info!("Using array position from Cotter.");
            LatLngHeight {
                longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
                latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
                height_metres: COTTER_MWA_HEIGHT_METRES,
            }
        } else {
            info!("Using default MWA array position.");
            LatLngHeight::new_mwa()
        };
        prep_ctx.phase_centre = match (
            matches
                .values_of_t::<f64>("phase-centre")
                .map(|v| (v[0], v[1])),
            matches.is_present("pointing-centre"),
        ) {
            (Err(err), _) if err.kind() != ArgumentNotFound => return Err(err.into()),
            (Ok(_), true) => {
                unreachable!("--phase-centre conflicts with --pointing-centre, enforced by clap");
            }
            (Ok((ra, dec)), _) => RADec::new(ra.to_radians(), dec.to_radians()),
            (_, true) => RADec::from_mwalib_tile_pointing(&corr_ctx.metafits_context),
            _ => RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context),
        };
        prep_ctx.correct_cable_lengths = {
            let no_cable_delays = matches.is_present("no-cable-delay");
            let cable_delays_applied = corr_ctx.metafits_context.cable_delays_applied;
            debug!(
                "cable corrections: applied={}, desired={}",
                cable_delays_applied, !no_cable_delays
            );
            !cable_delays_applied && !no_cable_delays
        };
        prep_ctx.correct_digital_gains = !matches.is_present("no-digital-gains");
        prep_ctx.passband_gains = match matches.value_of("passband-gains") {
            None | Some("none") => None,
            Some("jake") => Some(PFB_JAKE_2022_200HZ.to_vec()),
            Some("cotter") => Some(PFB_COTTER_2014_10KHZ.to_vec()),
            Some(option) => panic!("unknown option for --passband-gains: {}", option),
        };
        prep_ctx.correct_geometry = {
            let no_geometric_delays = matches.is_present("no-geometric-delay");
            let geometric_delays_applied = corr_ctx.metafits_context.geometric_delays_applied;
            debug!(
                "geometric corrections: applied={:?}, desired={}",
                geometric_delays_applied, !no_geometric_delays
            );
            matches!(geometric_delays_applied, GeometricDelaysApplied::No) && !no_geometric_delays
        };
        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                prep_ctx.aoflagger_strategy = if matches.is_present("no-rfi") {
                    None
                } else {
                    match matches.value_of_t("aoflagger-strategy") {
                        Err(err) if err.kind() != ArgumentNotFound => return Err(err.into()),
                        Ok(strategy) => Some(strategy),
                        Err(_) => Some(unsafe {
                            cxx_aoflagger_new().FindStrategyFileMWA()
                        }),
                    }
                };
            }
        }
        Ok(prep_ctx)
    }

    /// Parse an iterator of arguments, `args` into a `BirliContext`.
    ///
    /// # Errors
    ///
    /// Can raise:
    /// - `clap::Error` if clap cannot parse `args`
    /// - `mwalib::MwalibError` if mwalib can't open the input files.
    /// - `BirliError::CLIError` if the arguments are invalid.
    pub fn from_args<I, T>(args: I) -> Result<Self, BirliError>
    where
        I: IntoIterator<Item = T> + Debug,
        T: Into<OsString> + Clone,
    {
        debug!("args:\n{:?}", &args);

        let matches = Self::get_matches(args)?;
        trace!("arg matches:\n{:?}", &matches);

        for unimplemented_option in &[
            "flag-init",
            "flag-init-steps",
            "flag-end",
            "flag-end-steps",
            "flag-edge-width",
            "flag-edge-chans",
            "flag-dc",
            "no-flag-dc",
            "no-sel-autos",
            "no-sel-flagged-ants",
            "sel-ants",
        ] {
            assert!(
                !matches.is_present(unimplemented_option),
                "option not yet implemented: --{}",
                unimplemented_option
            );
        }

        for untested_option in &[
            "flag-times",
            "flag-coarse-chans",
            "flag-fine-chans",
            "flag-autos",
            "no-flag-metafits",
            "flag-antennas",
            "time-chunk",
            "max-memory",
        ] {
            if matches.is_present(untested_option) {
                warn!(
                    "option does not have full test coverage, use with caution: --{}",
                    untested_option
                );
            }
        }

        let io_ctx = Self::parse_io_matches(&matches);
        let corr_ctx = io_ctx.get_corr_ctx()?;
        debug!("mwalib correlator context:\n{}", &corr_ctx);
        let vis_sel = Self::parse_vis_sel_matches(&corr_ctx, &matches)?;
        let flag_ctx = Self::parse_flag_matches(&corr_ctx, &matches)?;
        let prep_ctx = Self::parse_prep_matches(&matches, &corr_ctx)?;
        let (avg_time, avg_freq) = Self::parse_avg_matches(&matches, &corr_ctx)?;
        let num_timesteps_per_chunk =
            Self::parse_chunk_matches(&corr_ctx, &matches, avg_time, &vis_sel)?;

        let result = Self {
            corr_ctx,
            prep_ctx,
            vis_sel,
            flag_ctx,
            io_ctx,
            avg_time,
            avg_freq,
            num_timesteps_per_chunk,
        };

        info!("{}", &result);

        if matches.is_present("dry-run") {
            return Err(DryRun {});
        }

        Ok(result)
    }

    /// Read, Preprocess and write corrected visibilities chunks.
    ///
    /// # Errors
    ///
    /// can raise:
    /// - `BirliError::BadArrayShape` if the shape of the calibartion solutions
    ///     is incompatible with the visibility shape.
    /// - preprocessing errors
    pub fn run(self) -> Result<HashMap<String, Duration>, BirliError> {
        let BirliContext {
            corr_ctx,
            mut prep_ctx,
            vis_sel,
            flag_ctx,
            io_ctx,
            avg_time,
            avg_freq,
            num_timesteps_per_chunk,
        } = self;

        // used to time large operations
        let mut durations = HashMap::<String, Duration>::new();

        // ////////// //
        // Prepare IO //
        // ////////// //

        let vis_ctx = VisContext::from_mwalib(
            &corr_ctx,
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            &vis_sel.baseline_idxs,
            avg_time,
            avg_freq,
        );

        // TODO: move phase_centre, array_pos out of prep_ctx
        let obs_ctx = ObsContext {
            phase_centre: prep_ctx.phase_centre,
            array_pos: prep_ctx.array_pos,
            ..ObsContext::from_mwalib(&corr_ctx.metafits_context)
        };

        let mwa_ctx = MwaObsContext::from_mwalib(&corr_ctx.metafits_context);

        prep_ctx.calsols = if let Some(ref calsol_file) = io_ctx.aocalsols_in {
            let calsols = AOCalSols::read_andre_binary(calsol_file).unwrap();
            assert!(
                calsols.di_jones.dim().0 == 1,
                "only 1 timeblock must be supplied for calsols. \
                Instead found {} timeblocks. dimensions {:?}",
                calsols.di_jones.dim().1,
                calsols.di_jones.dim()
            );
            let calsol_chans = calsols.di_jones.dim().2;
            if calsol_chans % corr_ctx.num_coarse_chans != 0 {
                return Err(BirliError::BadArrayShape {
                    argument: format!("io_ctx.aocalsols_in={}", calsol_file),
                    expected: format!(
                        "a multiple of metafits_num_coarse_chans={}",
                        corr_ctx.metafits_context.num_metafits_coarse_chans
                    ),
                    received: format!("{}", calsol_chans),
                    function: "BirliContext::run".into(),
                });
            }
            let num_calsol_fine_chans_per_coarse = calsol_chans / corr_ctx.num_coarse_chans;
            Some(
                calsols
                    .di_jones
                    .index_axis(Axis(0), 0)
                    .slice(s![
                        ..,
                        (vis_sel.coarse_chan_range.start * num_calsol_fine_chans_per_coarse)
                            ..(vis_sel.coarse_chan_range.end * num_calsol_fine_chans_per_coarse)
                    ])
                    .to_owned(),
            )
        } else {
            None
        };

        let mut uvfits_writer = io_ctx.uvfits_out.map(|uvfits_out| {
            with_increment_duration!(durations, "init", {
                UvfitsWriter::from_marlu(
                    uvfits_out,
                    &vis_ctx,
                    Some(obs_ctx.array_pos),
                    obs_ctx.phase_centre,
                    obs_ctx.name.clone(),
                )
                .expect("unable to initialize uvfits writer")
            })
        });
        let mut ms_writer = io_ctx.ms_out.map(|ms_out| {
            let writer =
                MeasurementSetWriter::new(ms_out, obs_ctx.phase_centre, Some(obs_ctx.array_pos));
            with_increment_duration!(durations, "init", {
                writer
                    .initialize_mwa(&vis_ctx, &obs_ctx, &mwa_ctx, &vis_sel.coarse_chan_range)
                    .expect("unable to initialize ms writer");
            });
            writer
        });

        let gpubox_ids = corr_ctx.coarse_chans[vis_sel.coarse_chan_range.clone()]
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect::<Vec<_>>();
        let mut flag_file_set = io_ctx.flag_template.map(|flag_template| {
            FlagFileSet::new(&flag_template, &gpubox_ids, corr_ctx.mwa_version)
                .expect("cannot create flag file writer")
        });

        // //////// //
        // Chunking //
        // //////// //

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let chunk_size = if let Some(steps) = num_timesteps_per_chunk {
            steps
        } else {
            vis_sel.timestep_range.len()
        };

        // Allocate our big arrays once, reuse them for each chunk unless the chunk shape changes
        let chunk_vis_sel = VisSelection {
            timestep_range: (vis_sel.timestep_range.start
                ..vis_sel.timestep_range.start + chunk_size),
            ..vis_sel.clone()
        };
        let mut jones_array = chunk_vis_sel.allocate_jones(fine_chans_per_coarse)?;
        let mut flag_array = chunk_vis_sel.allocate_flags(fine_chans_per_coarse)?;
        let mut weight_array = chunk_vis_sel.allocate_weights(fine_chans_per_coarse)?;

        for mut timestep_chunk in &vis_sel.timestep_range.clone().chunks(chunk_size) {
            let chunk_first_timestep = timestep_chunk.next().expect("zero-sized chunk");
            let chunk_vis_sel = VisSelection {
                timestep_range: (chunk_first_timestep
                    ..(timestep_chunk.last().unwrap_or(chunk_first_timestep) + 1)),
                ..vis_sel.clone()
            };
            if num_timesteps_per_chunk.is_some() {
                info!(
                    "processing timestep chunk {:?} of {:?} % {}",
                    chunk_vis_sel.timestep_range,
                    vis_sel.timestep_range.clone(),
                    chunk_size
                );
            }

            // only reallocate arrays if the chunk dimensions have changed.
            let chunk_dims = chunk_vis_sel.get_shape(fine_chans_per_coarse);
            if jones_array.dim() != chunk_dims {
                jones_array = chunk_vis_sel.allocate_jones(fine_chans_per_coarse)?;
            }
            if flag_array.dim() != chunk_dims {
                flag_array = chunk_vis_sel.allocate_flags(fine_chans_per_coarse)?;
            }
            if weight_array.dim() != chunk_dims {
                weight_array = chunk_vis_sel.allocate_weights(fine_chans_per_coarse)?;
            }

            // populate flags
            flag_ctx.set_flags(
                &mut flag_array,
                &chunk_vis_sel.timestep_range,
                &chunk_vis_sel.coarse_chan_range,
                &chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )?;

            // populate visibilities
            with_increment_duration!(
                durations,
                "read",
                chunk_vis_sel.read_mwalib(
                    &corr_ctx,
                    jones_array.view_mut(),
                    flag_array.view_mut(),
                    prep_ctx.draw_progress,
                )?
            );

            // populate weights
            weight_array.fill(vis_ctx.weight_factor() as f32);

            prep_ctx.preprocess(
                &corr_ctx,
                &mut jones_array,
                &mut weight_array,
                &mut flag_array,
                &mut durations,
                &chunk_vis_sel,
            )?;

            // output flags (before averaging)
            if let Some(flag_file_set) = flag_file_set.as_mut() {
                with_increment_duration!(
                    durations,
                    "write",
                    flag_file_set
                        .write_flag_array(&corr_ctx, &flag_array, &gpubox_ids)
                        .expect("unable to write flags")
                );
            }

            // bake flags into weights
            for (weight, flag) in izip!(weight_array.iter_mut(), flag_array.iter()) {
                *weight = if *flag {
                    -(*weight).abs()
                } else {
                    (*weight).abs()
                } as f32;
            }

            let chunk_vis_ctx = VisContext::from_mwalib(
                &corr_ctx,
                &chunk_vis_sel.timestep_range,
                &chunk_vis_sel.coarse_chan_range,
                &chunk_vis_sel.baseline_idxs,
                avg_time,
                avg_freq,
            );

            let ant_positions_geodetic: Vec<_> = obs_ctx.ant_positions_geodetic().collect();

            // output uvfits
            if let Some(uvfits_writer) = uvfits_writer.as_mut() {
                with_increment_duration!(
                    durations,
                    "write",
                    uvfits_writer
                        .write_vis_marlu(
                            jones_array.view(),
                            weight_array.view(),
                            &chunk_vis_ctx,
                            &ant_positions_geodetic,
                            prep_ctx.draw_progress,
                        )
                        .expect("unable to write uvfits")
                );
            }

            // output ms
            if let Some(ms_writer) = ms_writer.as_mut() {
                with_increment_duration!(
                    durations,
                    "write",
                    ms_writer
                        .write_vis_marlu(
                            jones_array.view(),
                            weight_array.view(),
                            &chunk_vis_ctx,
                            &ant_positions_geodetic,
                            prep_ctx.draw_progress,
                        )
                        .expect("unable to write ms")
                );
            }
        }

        // Finalise the uvfits writer.
        if let Some(uvfits_writer) = uvfits_writer {
            with_increment_duration!(
                durations,
                "write",
                uvfits_writer
                    .write_ants_from_mwalib(&corr_ctx.metafits_context)
                    .expect("couldn't write antenna table to uvfits")
            );
        };

        Ok(durations)
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
mod tests {
    use tempfile::tempdir;

    use crate::{test_common::get_1254670392_avg_paths, BirliContext};

    #[test]
    fn test_birli_context_display_doesnt_crash() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("high_2019B_2458765_EOR0"));
        assert!(display.contains("Will not correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will flag with aoflagger"));
        assert!(display.contains("Will not correct geometry"));
    }
}

#[cfg(test)]
mod argparse_tests {
    use crate::{error::BirliError, test_common::get_1254670392_avg_paths, BirliContext};

    #[test]
    fn test_parse_missing_input() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // no gpufits
        let args = vec!["birli", "-m", metafits_path];

        match BirliContext::from_args(&args) {
            Err(BirliError::ClapError(inner)) => assert!(matches!(
                inner.kind(),
                clap::error::ErrorKind::MissingRequiredArgument { .. }
            )),
            Err(e) => panic!("expected missing required argument error, not {}", e),
            Ok(_) => panic!("expected error, but got Ok(_)"),
        }

        // no metafits
        let mut args = vec!["birli"];
        args.extend_from_slice(&gpufits_paths);

        match BirliContext::from_args(&args) {
            Err(BirliError::ClapError(inner)) => assert!(matches!(
                inner.kind(),
                clap::error::ErrorKind::MissingRequiredArgument { .. }
            )),
            Err(e) => panic!("expected missing required argument error, not {}", e),
            Ok(_) => panic!("expected error, but got Ok(_)"),
        }
    }

    #[test]
    fn test_parse_invalid_input() {
        let (metafits_path, mut gpufits_paths) = get_1254670392_avg_paths();

        // test nonexistent metafits

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", "nonexistent.metafits"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::MwalibError(_))
        ));

        gpufits_paths[0] = "nonexistent.fits";

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::MwalibError(_))
        ));
    }

    #[test]
    fn test_parse_invalid_time_selection() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "0", "999999"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "2", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_selection() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "1", "2"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { vis_sel, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(vis_sel.timestep_range, 1..3);
    }

    #[test]
    fn test_parse_invalid_time_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-times", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-times", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.timestep_flags[2]);
        assert!(!flag_ctx.timestep_flags[1]);
    }

    #[test]
    fn test_parse_invalid_coarse_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-coarse-chans", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_coarse_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-coarse-chans", "2",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.coarse_chan_flags[2]);
        assert!(!flag_ctx.coarse_chan_flags[1]);
    }

    #[test]
    fn test_parse_invalid_fine_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-fine-chans", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_fine_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-fine-chans", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.fine_chan_flags[2]);
        assert!(!flag_ctx.fine_chan_flags[1]);
    }

    #[test]
    fn test_parse_invalid_antenna_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--flag-antennas", "0", "9999",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_antenna_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-antennas", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.antenna_flags[2]);
        assert!(!flag_ctx.antenna_flags[1]);
    }

    #[test]
    fn test_parse_no_flag_metafits_and_valid_antenna_flag() {
        let (_, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            // in this metafits file, all ants are flagged.
            "-m", "tests/data/1254670392_avg/1254670392.metafits",
            "--no-flag-metafits",
            "--flag-antennas", "2",
            "--",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.antenna_flags[2]);
        assert!(!flag_ctx.antenna_flags[1]);
    }

    #[test]
    fn test_parse_flag_autos() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--flag-autos"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.autos);
    }

    #[test]
    fn test_parse_invalid_avg_time() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-factor", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-res", "0.5"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-time-res", "4s"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::ClapError(_))
        ));
    }

    #[test]
    fn test_parse_invalid_avg_freq() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-factor", "0"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "0.5"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "4s"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::ClapError(_))
        ));
    }

    #[test]
    fn test_parse_valid_avg_freq() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--avg-freq-res", "120"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { avg_freq, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(avg_freq, 3);
    }

    #[test]
    fn test_parse_invalid_time_chunk() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // test when time_chunk is not a multiple of avg_time
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--avg-time-factor", "2",
            "--time-chunk", "1",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_time_chunk() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec!["birli", "-m", metafits_path, "--time-chunk", "2"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            num_timesteps_per_chunk,
            ..
        } = BirliContext::from_args(&args).unwrap();

        assert_eq!(num_timesteps_per_chunk, Some(2));
    }

    #[test]
    fn test_parse_invalid_max_memory() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        // test when max mem is less than 1 byte
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory",
            "0.0000000000000001",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

        // test when max mem can't fit a single averaged timestep
        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory", "0.32",
            "--avg-time-factor", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));
    }

    #[test]
    fn test_parse_valid_max_memory() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--max-memory", "0.64",
            "--sel-time", "0", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext {
            num_timesteps_per_chunk,
            ..
        } = BirliContext::from_args(&args).unwrap();

        assert_eq!(num_timesteps_per_chunk, Some(2));
    }
}

/// These are basically integration tests, but if I put them in a separate
/// module, then I don't get the coverage. :(
/// All these tests require aoflagger because they either require flagging or
/// use the --no-rfi option.
/// TODO: get unit test coverage to the point where this can be moved to a unit
/// test module.
#[cfg(test)]
#[cfg(feature = "aoflagger")]
mod tests_aoflagger {
    use std::path::PathBuf;

    use float_cmp::F32Margin;
    use marlu::{
        rubbl_casatables::{Table, TableOpenMode},
        RADec,
    };
    use tempfile::tempdir;

    use crate::{
        test_common::{compare_ms_with_csv, compare_uvfits_with_csv, get_1254670392_avg_paths},
        BirliContext,
    };

    #[test]
    fn compare_cotter_uvfits_nocorrect_rfi() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default(),
            false,
            false,
        );
    }

    #[test]
    fn compare_cotter_uvfits_nocorrect_norfi_timechunk1() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.chunked.uvfits");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--time-chunk", "1",
            "--sel-time", "0", "2"
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert_eq!(birli_ctx.prep_ctx.aoflagger_strategy, None);
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(1));

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    #[test]
    fn compare_cotter_uvfits_geom_cable_rfi() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.uvfits.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
            false,
        );
    }

    /// Test corrections using arbitrary phase centre.
    /// data generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -centre 00h00m00.0s 00d00m00.0s \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// ```bash
    /// python tests/data/dump_uvfits.py \
    ///     tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits \
    ///     --timestep-limit=12 --baseline-limit=12 --dump-mode=vis-only \
    ///     --dump-csv=tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits.csv
    /// ```
    #[test]
    fn compare_cotter_uvfits_geom_cable_rfi_phase_custom() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase0.uvfits.csv",
        );

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--phase-centre", "0.0", "0.0",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.prep_ctx.phase_centre, RADec { ra: 0., dec: 0. });
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            &uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
            false,
            false,
        );
    }

    /// Test corrections using pointing centre as phase centre.
    /// data generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -usepcentre \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// Then the following CASA commands
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_geom_cable_rfi_phase_pointing() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.corrected.phase-point.ms.csv",
        );

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--pointing-centre",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--sel-time", "0", "1",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        let pointing_centre =
            RADec::from_mwalib_tile_pointing(&birli_ctx.corr_ctx.metafits_context);
        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.prep_ctx.phase_centre, pointing_centre);
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
            true,
            true,
        );
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms \
    ///  -allowmissing \
    ///  -edgewidth 0 \
    ///  -endflag 0 \
    ///  -initflag 0 \
    ///  -noantennapruning \
    ///  -nocablelength \
    ///  -norfi \
    ///  -nogeom \
    ///  -noflagautos \
    ///  -noflagdcchannels \
    ///  -nosbgains \
    ///  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///  -nostats \
    ///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///  -full-apply tests/data/1254670392_avg/1254690096.bin \
    ///  tests/data/1254670392_avg/*gpubox*.fits
    /// ```
    ///
    /// then casa
    ///
    /// ```bash
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_nocorrect_norfi_cal() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.cal.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.ms.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal", "tests/data/1254670392_avg/1254690096.bin",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    #[test]
    /// Handle when calibration solution is provided with 24 channels, but a subset of channels are provided
    fn compare_cotter_ms_none_norfi_cal_partial() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.cal.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.norfi.cal.partial.ms.csv",
        );

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal", "tests/data/1254670392_avg/1254690096.bin",
        ];
        args.extend_from_slice(&gpufits_paths[21..]);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default(),
            true,
            false,
        );
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms \
    ///  -allowmissing \
    ///  -edgewidth 0 \
    ///  -endflag 0 \
    ///  -initflag 0 \
    ///  -noantennapruning \
    ///  -nocablelength \
    ///  -norfi \
    ///  -nogeom \
    ///  -noflagautos \
    ///  -noflagdcchannels \
    ///  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///  -nostats \
    ///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///  tests/data/1254670392_avg/*gpubox*.fits
    /// ```
    ///
    /// then casa
    ///
    /// ```bash
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_norfi_nopfb() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.nopfb.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.norfi.nopfb.ms.csv");

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct digital gains"));
        assert!(display.contains("Will not flag with aoflagger"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(7e-5),
            false,
            true,
        );
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms \
    ///  -allowmissing \
    ///  -edgewidth 0 \
    ///  -endflag 0 \
    ///  -initflag 0 \
    ///  -noantennapruning \
    ///  -nocablelength \
    ///  -norfi \
    ///  -nogeom \
    ///  -noflagautos \
    ///  -noflagdcchannels \
    ///  -nosbgains \
    ///  -sbpassband tests/data/subband-passband-32ch-cotter.txt \
    ///  -nostats \
    ///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///  tests/data/1254670392_avg/*gpubox*.fits
    /// ```
    ///
    /// then casa
    ///
    /// ```bash
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_norfi_nodigital_pfb_cotter_40() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.nodigital.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.pfb-cotter-40.ms.csv",
        );

        #[rustfmt::skip]
        let args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--no-digital-gains",
            "--pfb-gains", "cotter",
            "--emulate-cotter",
            gpufits_paths[23],
            gpufits_paths[22],
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(matches!(birli_ctx.prep_ctx.passband_gains, Some(_)));
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct coarse pfb passband gains"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-2),
            false,
            true,
        );
    }

    /// Test generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.corrected.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// then the following casa commands:
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.corrected.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_corrected() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.corrected.ms.csv");

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("Will correct cable lengths"));
        assert!(display.contains("Will not correct digital gains"));
        assert!(display.contains("Will not correct coarse pfb passband gains"));
        assert!(display.contains("Will flag with aoflagger"));
        assert!(display.contains("Will correct geometry"));

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
            false,
        );
    }

    /// Test generated with:
    ///
    /// ```bash
    /// cotter \
    ///   -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///   -o tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms \
    ///   -allowmissing \
    ///   -edgewidth 0 \
    ///   -endflag 0 \
    ///   -initflag 0 \
    ///   -noantennapruning \
    ///   -nocablelength \
    ///   -nogeom \
    ///   -noflagautos \
    ///   -noflagdcchannels \
    ///   -nosbgains \
    ///   -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///   -nostats \
    ///   -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///   -timeres 4 \
    ///   -freqres 160 \
    ///   tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```
    ///
    /// then the following casa commands:
    ///
    /// ```python
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms/')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-res", "4",
            "--avg-freq-res", "160",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }

    /// Same as above but with factors instead of resolution
    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz_factors() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }

    #[test]
    fn compare_cotter_ms_none_avg_4s_160khz_max_mem() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.ms.csv");

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-M", ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--sel-time", "0", "2",
            "--time-chunk", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(2));
        assert_eq!(birli_ctx.vis_sel.timestep_range, 0..3);

        birli_ctx.run().unwrap();

        let main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        assert_eq!(main_table.n_rows(), 2 * 8256);

        compare_ms_with_csv(
            &ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
            false,
        );
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger_flagset {
    use crate::{io::mwaf::FlagFileSet, BirliContext};
    use itertools::izip;
    use marlu::mwalib::CorrelatorContext;
    use tempfile::tempdir;

    macro_rules! assert_flagsets_eq {
        ($context:expr, $left_flagset:expr, $right_flagset:expr, $gpubox_ids:expr) => {
            let num_baselines = $context.metafits_context.num_baselines;
            let num_flags_per_row = $context.metafits_context.num_corr_fine_chans_per_coarse;
            let num_common_timesteps = $context.num_common_timesteps;
            let num_rows = num_common_timesteps * num_baselines;
            let num_flags_per_timestep = num_baselines * num_flags_per_row;

            assert!(num_baselines > 0);
            assert!(num_rows > 0);
            assert!(num_flags_per_row > 0);

            let right_chan_header_flags_raw =
            $right_flagset.read_chan_header_flags_raw().unwrap();

            let left_chan_header_flags_raw = $left_flagset.read_chan_header_flags_raw().unwrap();

            for gpubox_id in $gpubox_ids {
                let (left_header, left_flags) = left_chan_header_flags_raw.get(&gpubox_id).unwrap();
                let (right_header, right_flags) =
                    right_chan_header_flags_raw.get(&gpubox_id).unwrap();
                assert_eq!(left_header.obs_id, right_header.obs_id);
                assert_eq!(left_header.num_channels, right_header.num_channels);
                assert_eq!(left_header.num_ants, right_header.num_ants);
                // assert_eq!(left_header.num_common_timesteps, right_header.num_common_timesteps);
                assert_eq!(left_header.num_timesteps, num_common_timesteps);
                assert_eq!(left_header.num_pols, right_header.num_pols);
                assert_eq!(left_header.gpubox_id, right_header.gpubox_id);
                assert_eq!(left_header.bytes_per_row, right_header.bytes_per_row);
                // assert_eq!(left_header.num_rows, right_header.num_rows);
                assert_eq!(left_header.num_rows, num_rows);

                // assert_eq!(left_flags.len(), right_flags.len());
                assert_eq!(
                    left_flags.len(),
                    num_common_timesteps * num_baselines * num_flags_per_row
                );

                izip!(
                    left_flags.chunks(num_flags_per_timestep),
                    right_flags.chunks(num_flags_per_timestep)
                ).enumerate().for_each(|(common_timestep_idx, (left_timestep_chunk, right_timestep_chunk))| {
                    izip!(
                        $context.metafits_context.baselines.iter(),
                        left_timestep_chunk.chunks(num_flags_per_row),
                        right_timestep_chunk.chunks(num_flags_per_row)
                    ).enumerate().for_each(|(baseline_idx, (baseline, left_baseline_chunk, right_baseline_chunk))| {
                        if baseline.ant1_index == baseline.ant2_index {
                            return
                        }

                        assert_eq!(
                            left_baseline_chunk, right_baseline_chunk,
                            "flag chunks for common timestep {}, baseline {} (ants {}, {}) do not match! \nbirli:\n{:?}\ncotter:\n{:?}",
                            common_timestep_idx, baseline_idx, baseline.ant1_index, baseline.ant2_index, left_baseline_chunk, right_baseline_chunk
                        )
                    });
                });
            }
        };
    }

    #[test]
    fn aoflagger_outputs_flags() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "none",
            "-f", mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(matches!(birli_ctx.prep_ctx.passband_gains, None));
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.flag_template,
            Some(mwaf_path_template.to_str().unwrap().into())
        );

        let corr_ctx =
            CorrelatorContext::new(&birli_ctx.io_ctx.metafits_in, &birli_ctx.io_ctx.gpufits_in)
                .unwrap();

        birli_ctx.run().unwrap();

        let gpubox_ids: Vec<usize> = corr_ctx
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| corr_ctx.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let mut birli_flag_file_set = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();

        let mut cotter_flag_file_set = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();

        assert_flagsets_eq!(
            &corr_ctx,
            birli_flag_file_set,
            cotter_flag_file_set,
            gpubox_ids
        );
    }

    #[test]
    #[should_panic]
    fn aoflagger_outputs_flags_chunked() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--time-chunk", "1",
            "-f", mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        BirliContext::from_args(&args).unwrap();
    }
}
