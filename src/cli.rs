//! Command Line Interface helpers for Birli

use crate::{
    context_to_jones_array,
    error::{BirliError, BirliError::DryRun, CLIError::InvalidCommandLineArgument},
    flags::FlagContext,
    flags::{add_dimension, flag_to_weight_array, get_weight_factor},
    io::IOContext,
    io::{aocal::AOCalSols, WriteableVis},
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        hifitime::Epoch,
        mwalib::{CorrelatorContext, GeometricDelaysApplied},
        precession::{precess_time, PrecessionInfo},
        LatLngHeight, RADec,
    },
    marlu::{
        io::{ms::MeasurementSetWriter, VisWritable},
        ndarray::s,
    },
    passband_gains::{PFB_COTTER_2014_10KHZ, PFB_JAKE_2022_200HZ},
    with_increment_duration, Axis, Complex, FlagFileSet, PreprocessContext, UvfitsWriter,
    VisSelection,
};
use cfg_if::cfg_if;
use clap::{arg, command, PossibleValue, ValueHint::FilePath};
use itertools::Itertools;
use log::{debug, info, trace, warn};
use prettytable::{cell, format as prettyformat, row, table};
use std::{collections::HashMap, time::Duration};
use std::{
    env,
    ffi::OsString,
    fmt::{Debug, Display},
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
            if self.avg_time != 1 {
                format!(" ({}x)", self.avg_time)
            } else {
                "".into()
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
            if self.avg_freq != 1 {
                format!(" ({}x)", self.avg_freq)
            } else {
                "".into()
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

        let show_timestep_table = true;

        writeln!(
            f,
            "Timestep details (all={}, provided={}, common={}, good={}, select={}, flag={}):{}",
            self.corr_ctx.num_timesteps,
            self.corr_ctx.num_provided_timesteps,
            self.corr_ctx.num_common_timesteps,
            self.corr_ctx.num_common_good_timesteps,
            self.vis_sel.timestep_range.len(),
            timestep_flag_idxs.len(),
            if show_timestep_table {
                format!("\n{}", timestep_table)
            } else {
                "".into()
            }
        )?;
        if !show_timestep_table {
            writeln!(
                f,
                "-> provided:    {:?}",
                self.corr_ctx.provided_timestep_indices
            )?;
            writeln!(
                f,
                "-> common:      {:?}",
                self.corr_ctx.common_timestep_indices
            )?;
            writeln!(
                f,
                "-> common good: {:?}",
                self.corr_ctx.common_good_timestep_indices
            )?;
            writeln!(f, "-> selected:    {:?}", self.vis_sel.timestep_range)?;
        }

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

        let show_coarse_chan_table = true;

        writeln!(
            f,
            "Coarse channel details (metafits={}, provided={}, common={}, good={}, select={}, flag={}):{}",
            self.corr_ctx.num_coarse_chans,
            self.corr_ctx.num_provided_coarse_chans,
            self.corr_ctx.num_common_coarse_chans,
            self.corr_ctx.num_common_good_coarse_chans,
            self.vis_sel.coarse_chan_range.len(),
            coarse_chan_flag_idxs.len(),
            if show_coarse_chan_table { format!("\n{}", coarse_chan_table) } else { "".into() }
        )?;

        if !show_coarse_chan_table {
            writeln!(
                f,
                "-> provided:    {:?}",
                self.corr_ctx.provided_coarse_chan_indices
            )?;
            writeln!(
                f,
                "-> common:      {:?}",
                self.corr_ctx.common_coarse_chan_indices
            )?;
            writeln!(
                f,
                "-> common good: {:?}",
                self.corr_ctx.common_good_coarse_chan_indices
            )?;
            writeln!(f, "-> selected:    {:?}", self.vis_sel.coarse_chan_range)?;
        }

        let mut ant_table = table!([
            "",
            "tile",
            "name",
            "north [m]",
            "east [m]",
            "height [m]",
            "f"
        ]);
        ant_table.set_format(*prettyformat::consts::FORMAT_CLEAN);

        for (ant_idx, ant) in self.corr_ctx.metafits_context.antennas.iter().enumerate() {
            let flagged = *self.flag_ctx.antenna_flags.get(ant_idx).unwrap_or(&false);
            let row = row![r =>
                format!("ant{}:", ant_idx),
                ant.tile_id,
                ant.tile_name,
                format!("{:.3}", ant.north_m),
                format!("{:.3}", ant.east_m),
                format!("{:.3}", ant.height_m),
                if flagged {"f"} else {""}
            ];
            ant_table.add_row(row);
        }

        debug!(
            "Antenna details (all={}, flag={}):{}",
            self.corr_ctx.metafits_context.num_ants,
            self.flag_ctx
                .antenna_flags
                .iter()
                .enumerate()
                .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
                .count(),
            format!("\n{}", ant_table)
        );

        // let show_baseline_table = false;

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
        let mem_per_timestep_gib = (num_sel_chans
            * num_sel_baselines
            * num_sel_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())) as f64
            / 1024.0_f64.powi(3);

        writeln!(
            f,
            "Estimated memory usage per timestep =           {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_chans,
            num_sel_baselines,
            num_sel_pols,
            std::mem::size_of::<Complex<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib,
        )?;

        if let Some(num_timesteps) = self.num_timesteps_per_chunk {
            writeln!(
                f,
                "Estimated memory per chunk          = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
                num_timesteps,
                num_sel_chans,
                num_sel_baselines,
                num_sel_pols,
                std::mem::size_of::<Complex<f32>>(),
                std::mem::size_of::<f32>(),
                std::mem::size_of::<bool>(),
                mem_per_timestep_gib * num_timesteps as f64,
            )?;
        }

        writeln!(
            f,
            "Estimated memory selected           = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_sel_timesteps,
            num_sel_chans,
            num_sel_baselines,
            num_sel_pols,
            std::mem::size_of::<Complex<f32>>(),
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
    /// Parse an iterator of arguments, `args` into a `BirliContext`.
    ///
    /// # Errors
    ///
    /// TODO: parsing error handling
    ///
    /// Can raise:
    /// - `clap::Error` if clap cannot parse `args`
    /// - `mwalib::MwalibError` if mwalib can't open the input files.
    pub fn from_args<I, T>(args: I) -> Result<BirliContext, BirliError>
    where
        I: IntoIterator<Item = T>,
        T: Into<OsString> + Clone,
        I: Debug,
    {
        debug!("args:\n{:?}", &args);

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
                    .required(false),
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
        trace!("arg matches:\n{:?}", &matches);

        // // //
        // IO //
        // // //

        let io_ctx = IOContext {
            metafits_in: matches
                .value_of("metafits")
                .expect("--metafits is a required argument, enforced by clap")
                .into(),
            gpufits_in: matches
                .values_of("fits_paths")
                .expect("<PATHS> is a required argument, enforced by clap")
                .map(|s| s.into())
                .collect(),
            aocalsols_in: matches.value_of("apply-di-cal").map(|s| s.into()),
            uvfits_out: matches.value_of("uvfits-out").map(|s| s.into()),
            ms_out: matches.value_of("ms-out").map(|s| s.into()),
            flag_template: matches.value_of("flag-template").map(|s| s.into()),
        };

        let corr_ctx = io_ctx.get_corr_ctx()?;
        debug!("mwalib correlator context:\n{}", &corr_ctx);

        // ////////// //
        // Selections //
        // ////////// //

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

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
                clap::ErrorKind::ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        }

        // TODO: sel-ants, no-sel-flagged-ants, no-sel-autos

        // //////// //
        // Flagging //
        // //////// //

        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);

        // Timesteps

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
                clap::ErrorKind::ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };

        // TODO: flag-init, flag-steps, flag-end, flag-end-steps,

        // coarse channels
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
                clap::ErrorKind::ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };

        // fine channels
        // TODO: flag-edge-width, flag-edge-chans, flag-dc, no-flag-dc,
        match matches.values_of_t::<usize>("flag-fine-chans") {
            Ok(fine_chan_idxs) => {
                for (value_idx, &fine_chan_idx) in fine_chan_idxs.iter().enumerate() {
                    if fine_chan_idx >= corr_ctx.metafits_context.num_corr_fine_chans_per_coarse {
                        return Err(BirliError::CLIError(InvalidCommandLineArgument {
                            option: "--flag-fine-chans <CHANS>...".into(),
                            expected: format!(
                                "fine_chan_idx < num_fine_chans={}",
                                corr_ctx.metafits_context.num_corr_fine_chans_per_coarse
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
                clap::ErrorKind::ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };

        // Antennas
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
                clap::ErrorKind::ArgumentNotFound { .. } => {}
                _ => return Err(err.into()),
            },
        };

        // Baselines
        if matches.is_present("flag-autos") {
            flag_ctx.autos = true;
        }

        // ///////// //
        // Averaging //
        // ///////// //

        let int_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1e3;

        let avg_time: usize = match (
            matches.value_of("avg-time-factor"),
            matches.value_of("avg-time-res"),
        ) {
            (Some(_), Some(_)) => {
                panic!("you can't use --avg-time-factor and --avg-time-res at the same time");
            }
            (Some(factor_str), None) => factor_str.parse().unwrap_or_else(|_| {
                panic!(
                    "unable to parse --avg-time-factor \"{}\" as an unsigned integer",
                    factor_str
                )
            }),
            (_, Some(res_str)) => {
                let res = res_str.parse::<f64>().unwrap_or_else(|_| {
                    panic!("unable to parse --avg-time-res \"{}\" as a float", res_str)
                });
                let ratio = res / int_time_s;
                assert!(
                    ratio.is_finite() && ratio >= 1.0 && ratio.fract() < 1e-6,
                    "--avg-time-res {} must be an integer multiple of the input resolution, {}",
                    res,
                    int_time_s
                );
                ratio.round() as _
            }
            _ => 1,
        };

        let fine_chan_width_khz = corr_ctx.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;

        let avg_freq: usize = match (
            matches.value_of("avg-freq-factor"),
            matches.value_of("avg-freq-res"),
        ) {
            (Some(_), Some(_)) => {
                panic!("you can't use --avg-freq-factor and --avg-freq-res at the same time");
            }
            (Some(factor_str), None) => factor_str.parse().unwrap_or_else(|_| {
                panic!(
                    "unable to parse --avg-freq-factor \"{}\" as an unsigned integer",
                    factor_str
                )
            }),
            (_, Some(res_str)) => {
                let res = res_str.parse::<f64>().unwrap_or_else(|_| {
                    panic!("unable to parse --avg-freq-res \"{}\" as a float", res_str)
                });
                let ratio = res / fine_chan_width_khz;
                assert!(
                    ratio.is_finite() && ratio >= 1.0 && ratio.fract() < 1e-6,
                    "--avg-freq-res {} must be an integer multiple of the input resolution, {}",
                    res,
                    fine_chan_width_khz
                );
                ratio.round() as _
            }
            _ => 1,
        };

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let num_sel_timesteps = vis_sel.timestep_range.len();
        let num_sel_chans = vis_sel.coarse_chan_range.len() * fine_chans_per_coarse;
        let num_sel_baselines = vis_sel.baseline_idxs.len();
        let num_sel_pols = corr_ctx.metafits_context.num_visibility_pols;
        let bytes_per_timestep = num_sel_chans
            * num_sel_baselines
            * num_sel_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>());

        let num_timesteps_per_chunk: Option<usize> = match (
            matches.value_of("time-chunk"),
            matches.value_of("max-memory"),
        ) {
            (Some(_), Some(_)) => {
                // TODO: custom error type
                panic!("you can't use --time-chunk and --max-memory at the same time");
            }
            (Some(steps_str), None) => {
                // TODO: custom error type
                let steps = steps_str.parse().unwrap_or_else(|_| {
                    panic!(
                        "unable to parse --time-chunk \"{}\" as an unsigned integer",
                        steps_str
                    )
                });
                if steps % avg_time != 0 {
                    panic!(
                        "--time-chunk {} must be an integer multiple of the averaging factor, {}",
                        steps, avg_time
                    );
                }
                Some(steps)
            }
            (_, Some(mem_str)) => {
                // TODO: custom error type
                let max_memory_bytes = mem_str.parse::<f64>().unwrap_or_else(|_| {
                    panic!("unable to parse --max-memory \"{}\" as a float", mem_str)
                }) * 1024.0_f64.powi(3);
                if max_memory_bytes < 1.0 {
                    panic!(
                        "--max-memory must be at least 1 Byte, not {}B",
                        max_memory_bytes
                    );
                }
                let bytes_per_avg_time = bytes_per_timestep * avg_time;
                let num_bytes_total = num_sel_timesteps * bytes_per_timestep;
                if max_memory_bytes < num_bytes_total as f64 {
                    if max_memory_bytes < bytes_per_avg_time as f64 {
                        panic!(
                            "--max-memory ({} GiB) too small to fit a single averaged timestep ({} * {:.02} = {:.02} GiB)",
                            max_memory_bytes as f64 / 1024.0_f64.powi(3), avg_time, bytes_per_timestep as f64 / 1024.0_f64.powi(3), bytes_per_avg_time as f64 / 1024.0_f64.powi(3)
                        );
                    }
                    Some((max_memory_bytes / bytes_per_avg_time as f64).floor() as usize * avg_time)
                } else {
                    None
                }
            }
            _ => None,
        };

        // validate chunk size
        if let Some(chunk_size) = num_timesteps_per_chunk {
            if matches.value_of("flag-template").is_some() {
                panic!("chunking is not supported when writing .mwaf files using --flag-template");
            }
            info!("chunking output to {} timesteps per chunk", chunk_size);
        }

        // ////////////////// //
        // Correction Options //
        // ////////////////// //

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
            matches.values_of("phase-centre"),
            matches.is_present("pointing-centre"),
        ) {
            (Some(_), true) => {
                // TODO: custom error type
                panic!("--phase-centre can't be used with --pointing-centre");
            }
            (Some(mut values), _) => {
                // TODO: custom error type
                if let (Some(ra), Some(dec)) = (values.next(), values.next()) {
                    let ra = ra
                        .parse::<f64>()
                        .unwrap_or_else(|_| panic!("unable to parse RA {}", ra));
                    let dec = dec
                        .parse::<f64>()
                        .unwrap_or_else(|_| panic!("unable to parse DEC {}", dec));
                    debug!(
                        "Using phase centre from command line: RA={}, DEC={}",
                        ra, dec
                    );
                    RADec::new(ra.to_radians(), dec.to_radians())
                } else {
                    panic!("Unable to parse RADec. from --phase-centre");
                }
            }
            (_, true) => RADec::from_mwalib_tile_pointing(&corr_ctx.metafits_context),
            _ => RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context),
        };

        // cable delay corrections are enabled by default if they haven't aleady beeen applied.
        let no_cable_delays = matches.is_present("no-cable-delay");
        let cable_delays_applied = corr_ctx.metafits_context.cable_delays_applied;
        debug!(
            "cable corrections: applied={}, desired={}",
            cable_delays_applied, !no_cable_delays
        );
        prep_ctx.correct_cable_lengths = !cable_delays_applied && !no_cable_delays;

        // coarse channel digital gain corrections are enabled by default
        prep_ctx.correct_digital_gains = !matches.is_present("no-digital-gains");

        // coarse pfb passband corrections are enabled by default
        prep_ctx.passband_gains = match matches.value_of("passband-gains") {
            None | Some("none") => None,
            Some("jake") => Some(PFB_JAKE_2022_200HZ.to_vec()),
            Some("cotter") => Some(PFB_COTTER_2014_10KHZ.to_vec()),
            Some(option) => panic!("unknown option for --passband-gains: {}", option),
        };

        // geometric corrections are enabled by default if they haven't aleady beeen applied.
        let no_geometric_delays = matches.is_present("no-geometric-delay");
        let geometric_delays_applied = corr_ctx.metafits_context.geometric_delays_applied;
        debug!(
            "geometric corrections: applied={:?}, desired={}",
            geometric_delays_applied, !no_geometric_delays
        );
        prep_ctx.correct_geometry =
            matches!(geometric_delays_applied, GeometricDelaysApplied::No) && !no_geometric_delays;

        cfg_if! {
            if #[cfg(feature = "aoflagger")] {

                prep_ctx.aoflagger_strategy = if !matches.is_present("no-rfi") {
                    Some(matches.value_of("aoflagger-strategy").map(|s| s.into()).unwrap_or_else(|| unsafe {
                        cxx_aoflagger_new().FindStrategyFileMWA()
                    }))
                } else {
                    None
                };
            }
        }

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
            if matches.is_present(unimplemented_option) {
                panic!("option not yet implemented: --{}", unimplemented_option);
            }
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

        // ///////// //
        // Show info //
        // ///////// //

        let result = BirliContext {
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
    /// TODO: Errors docco
    pub fn run(self) -> Result<HashMap<String, Duration>, BirliError> {
        let BirliContext {
            corr_ctx,
            prep_ctx,
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

        let calsols_owned = io_ctx.aocalsols_in.map(|calsol_file| {
            let calsols = AOCalSols::read_andre_binary(calsol_file).unwrap();
            if calsols.di_jones.dim().0 != 1 {
                panic!("only 1 timeblock must be supplied for calsols. Instead found {} timeblocks. dimensions {:?}", calsols.di_jones.dim().1, calsols.di_jones.dim());
            }
            // calsols.di_jones.index_axis_move(Axis(0), 0)
            let calsol_chans = calsols.di_jones.dim().2;
            if calsol_chans % corr_ctx.num_coarse_chans != 0 {
                panic!(
                    "the number of calibration solution channels must be a multiple of the number of
                    coarse channels defined in the metafits {}. Instead found {}.
                    dimensions: {:?}",
                    corr_ctx.metafits_context.num_metafits_coarse_chans,
                    calsol_chans,
                    calsols.di_jones.dim());
            }
            let num_calsol_fine_chans_per_coarse = calsol_chans / corr_ctx.num_coarse_chans;
            calsols.di_jones
                .index_axis(Axis(0), 0)
                .slice(s![
                    ..,
                    (vis_sel.coarse_chan_range.start * num_calsol_fine_chans_per_coarse)
                    ..(vis_sel.coarse_chan_range.end * num_calsol_fine_chans_per_coarse)]
                ).to_owned()
        });

        let mut uvfits_writer = io_ctx.uvfits_out.map(|uvfits_out| {
            with_increment_duration!(durations, "init", {
                UvfitsWriter::from_mwalib(
                    uvfits_out,
                    &corr_ctx,
                    &vis_sel.timestep_range,
                    &vis_sel.coarse_chan_range,
                    &vis_sel.baseline_idxs,
                    Some(prep_ctx.array_pos),
                    Some(prep_ctx.phase_centre),
                    avg_time,
                    avg_freq,
                )
                .expect("couldn't initialise uvfits writer")
            })
        });
        let mut ms_writer = io_ctx.ms_out.map(|ms_out| {
            let writer =
                MeasurementSetWriter::new(ms_out, prep_ctx.phase_centre, Some(prep_ctx.array_pos));
            with_increment_duration!(durations, "init", {
                writer
                    .initialize_from_mwalib(
                        &corr_ctx,
                        &vis_sel.timestep_range,
                        &vis_sel.coarse_chan_range,
                        &vis_sel.baseline_idxs,
                        avg_time,
                        avg_freq,
                    )
                    .unwrap();
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

        let full_sel_timestep_range = vis_sel.timestep_range.clone();
        let chunk_size = if let Some(steps) = num_timesteps_per_chunk {
            steps
        } else {
            full_sel_timestep_range.len()
        };
        for mut timestep_chunk in &full_sel_timestep_range.clone().chunks(chunk_size) {
            let chunk_first_timestep = timestep_chunk.next().unwrap();
            let chunk_last_timestep = timestep_chunk.last().unwrap_or(chunk_first_timestep);
            let chunk_vis_sel = VisSelection {
                timestep_range: chunk_first_timestep..chunk_last_timestep + 1,
                ..vis_sel.clone()
            };
            if num_timesteps_per_chunk.is_some() {
                info!(
                    "processing timestep chunk {:?} of {:?} % {}",
                    chunk_vis_sel.timestep_range,
                    full_sel_timestep_range.clone(),
                    chunk_size
                );
            }
            let flag_array = flag_ctx.to_array(
                &chunk_vis_sel.timestep_range,
                &chunk_vis_sel.coarse_chan_range,
                chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            );

            let (mut jones_array, mut flag_array) = with_increment_duration!(
                durations,
                "read",
                context_to_jones_array(
                    &corr_ctx,
                    &chunk_vis_sel.timestep_range,
                    &chunk_vis_sel.coarse_chan_range,
                    Some(flag_array),
                    prep_ctx.draw_progress,
                )
                .unwrap()
            );

            // generate weights
            let weight_factor = get_weight_factor(&corr_ctx);
            let mut weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

            prep_ctx
                .preprocess(
                    &corr_ctx,
                    &mut jones_array,
                    &mut weight_array,
                    &mut flag_array,
                    &calsols_owned,
                    &mut durations,
                    &chunk_vis_sel,
                )
                .expect("unable to preprocess the chunk.");

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

            // let marlu_context = MarluVisContext::from_mwalib(
            //     &context,
            //     &chunk_timestep_range,
            //     &coarse_chan_range,
            //     &chunk_vis_sel.baseline_idxs,
            //     avg_time,
            //     avg_freq,
            // );

            // TODO: nothing actually uses the pol axis for flags and weights, so rip it out.
            let num_pols = corr_ctx.metafits_context.num_visibility_pols;
            let flag_array = add_dimension(flag_array.view(), num_pols);
            let weight_array = add_dimension(weight_array.view(), num_pols);

            // output uvfits
            if let Some(uvfits_writer) = uvfits_writer.as_mut() {
                with_increment_duration!(
                    durations,
                    "write",
                    uvfits_writer
                        .write_vis_mwalib(
                            jones_array.view(),
                            weight_array.view(),
                            flag_array.view(),
                            &corr_ctx,
                            &chunk_vis_sel.timestep_range,
                            &chunk_vis_sel.coarse_chan_range,
                            &chunk_vis_sel.baseline_idxs,
                            avg_time,
                            avg_freq,
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
                        .write_vis_mwalib(
                            jones_array.view(),
                            weight_array.view(),
                            flag_array.view(),
                            &corr_ctx,
                            &chunk_vis_sel.timestep_range,
                            &chunk_vis_sel.coarse_chan_range,
                            &chunk_vis_sel.baseline_idxs,
                            avg_time,
                            avg_freq,
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

    use crate::{test_common::*, BirliContext};

    #[test]
    fn test_birli_context_display_doesnt_crash() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );

        let display = format!("{}", &birli_ctx);
        assert!(display.contains("high_2019B_2458765_EOR0"));
        assert!(display.contains("Will not correct cable lengths"));
    }
}

#[cfg(test)]
mod argparse_tests {
    use crate::{error::BirliError, test_common::*, BirliContext};

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

        let mut args = vec!["birli", "-m", "nonexistent.metafits"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::MwalibError(_))
        ));

        gpufits_paths[0] = "nonexistent.fits";

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

        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "0", "999999"];
        args.extend_from_slice(&gpufits_paths);

        assert!(matches!(
            BirliContext::from_args(&args),
            Err(BirliError::CLIError(_))
        ));

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

        let mut args = vec!["birli", "-m", metafits_path, "--sel-time", "1", "2"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { vis_sel, .. } = BirliContext::from_args(&args).unwrap();

        assert_eq!(vis_sel.timestep_range, 1..3);
    }

    #[test]
    fn test_parse_invalid_time_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--flag-times",
            "0",
            "9999",
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

        let mut args = vec!["birli", "-m", metafits_path, "--flag-times", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.timestep_flags[2]);
        assert!(!flag_ctx.timestep_flags[1]);
    }

    #[test]
    fn test_parse_invalid_coarse_chan_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--flag-coarse-chans",
            "0",
            "9999",
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--flag-coarse-chans",
            "2",
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--flag-fine-chans",
            "0",
            "9999",
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

        let mut args = vec!["birli", "-m", metafits_path, "--flag-fine-chans", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.fine_chan_flags[2]);
        assert!(!flag_ctx.fine_chan_flags[1]);
    }

    #[test]
    fn test_parse_invalid_antenna_flag() {
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--flag-antennas",
            "0",
            "9999",
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

        let mut args = vec!["birli", "-m", metafits_path, "--flag-antennas", "2", "--"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.antenna_flags[2]);
        assert!(!flag_ctx.antenna_flags[1]);
    }

    #[test]
    fn test_parse_no_flag_metafits_and_valid_antenna_flag() {
        let (_, gpufits_paths) = get_1254670392_avg_paths();

        let mut args = vec![
            "birli",
            "-m",
            // in this metafits file, all ants are flagged.
            "tests/data/1254670392_avg/1254670392.metafits",
            "--no-flag-metafits",
            "--flag-antennas",
            "2",
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

        let mut args = vec!["birli", "-m", metafits_path, "--flag-autos"];
        args.extend_from_slice(&gpufits_paths);

        let BirliContext { flag_ctx, .. } = BirliContext::from_args(&args).unwrap();

        assert!(flag_ctx.autos);
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

    use crate::{test_common::*, BirliContext};

    #[test]
    fn compare_cotter_uvfits_nocorrect_rfi() {
        let tmp_dir = tempdir().unwrap();
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), false);
    }

    #[test]
    fn compare_cotter_uvfits_nocorrect_norfi_timechunk1() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.chunked.uvfits");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--time-chunk",
            "1",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(1));

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), true);
    }

    /// Test data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.cable.uvfits \
    ///  -allowmissing \
    ///  -edgewidth 0 \
    ///  -endflag 0 \
    ///  -initflag 0 \
    ///  -noantennapruning \
    ///  -nogeom \
    ///  -noflagautos \
    ///  -noflagdcchannels \
    ///  -nosbgains \
    ///  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
    ///  -nostats \
    ///  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
    ///  tests/data/1254670392_avg/1254670392*gpubox*.fits
    /// ```

    #[test]
    #[ignore = "slow"]
    fn compare_cotter_uvfits_cable_only_rfi() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.cable.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(birli_ctx.prep_ctx.correct_cable_lengths);
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
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
            false,
        );
    }

    #[test]
    // #[ignore = "slow"]
    fn compare_cotter_uvfits_geom_only_rfi() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.geom.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, Some(_)));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.num_timesteps_per_chunk, None);

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-digital-gains",
            "--phase-centre",
            "0.0",
            "0.0",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.uvfits_out,
            Some(uvfits_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--pointing-centre",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
            false,
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal",
            "tests/data/1254670392_avg/1254690096.bin",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(ms_path, expected_csv_path, F32Margin::default(), true);
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
            "--apply-di-cal",
            "tests/data/1254670392_avg/1254690096.bin",
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
        compare_ms_with_csv(ms_path, expected_csv_path, F32Margin::default(), true);
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert_eq!(birli_ctx.prep_ctx.passband_gains, None);
        assert!(birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(7e-5),
            false,
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--no-digital-gains",
            "--pfb-gains",
            "cotter",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        let birli_ctx = BirliContext::from_args(&args).unwrap();

        assert!(!birli_ctx.prep_ctx.correct_cable_lengths);
        assert!(matches!(birli_ctx.prep_ctx.passband_gains, Some(_)));
        assert!(!birli_ctx.prep_ctx.correct_digital_gains);
        assert!(!birli_ctx.prep_ctx.correct_geometry);
        assert!(matches!(birli_ctx.prep_ctx.aoflagger_strategy, None));
        assert_eq!(birli_ctx.io_ctx.metafits_in, metafits_path.to_string());
        assert_eq!(
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-2),
            false,
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-res",
            "4",
            "--avg-freq-res",
            "160",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor",
            "2",
            "--avg-freq-factor",
            "4",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);

        birli_ctx.run().unwrap();

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
        );
    }

    /// Same as above but forcing chunks by using a small --max-memory
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
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(2));
        // TODO: ???
        // assert_eq!( birli_ctx.vis_sel.timestep_range, 0..1 );

        birli_ctx.run().unwrap();

        let main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        assert_eq!(main_table.n_rows(), 2 * 8256);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
        );
    }

    /// test when time_chunk is not a multiple of avg_time
    #[test]
    #[should_panic]
    fn compare_cotter_ms_none_avg_4s_160khz_tiny_chunk() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.ms");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        env_logger::try_init().unwrap_or(());

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--sel-time", "0", "2",
            "--time-chunk", "1",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );
        assert_eq!(birli_ctx.avg_time, 2);
        assert_eq!(birli_ctx.avg_freq, 4);
        assert_eq!(birli_ctx.num_timesteps_per_chunk, Some(1));
        assert_eq!(birli_ctx.vis_sel.timestep_range, 0..1);

        birli_ctx.run().unwrap();
    }

    /// Data generated with
    ///
    /// ```bash
    /// cotter \
    ///  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
    ///  -o tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms \
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
    ///  -nostats \
    ///  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
    ///  tests/data/1254670392_avg/*gpubox*.fits
    /// ```
    ///
    /// then casa
    ///
    /// ```bash
    /// tb.open('tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms')
    /// exec(open('tests/data/casa_dump_ms.py').read())
    /// ```
    #[ignore = "Cotter doesn't correctly average passband gains"]
    #[test]
    fn compare_cotter_ms_none_norfi_nodigital_pfb_cotter() {
        let tmp_dir = tempdir().unwrap();
        let ms_path = tmp_dir.path().join("1254670392.none.norfi.nodigital.ms");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path = PathBuf::from(
            "tests/data/1254670392_avg/1254670392.cotter.none.norfi.nodigital.ms.csv",
        );

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-M",
            ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--no-rfi",
            "--no-digital-gains",
            "--pfb-gains",
            "cotter",
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
            birli_ctx.io_ctx.gpufits_in,
            gpufits_paths.map(|p| p.to_string())
        );
        assert_eq!(
            birli_ctx.io_ctx.ms_out,
            Some(ms_path.to_str().unwrap().into())
        );

        birli_ctx.run().unwrap();
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
            false,
        );
    }
}
