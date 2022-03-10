//! Command Line Interface helpers for Birli

use crate::{
    error::{BirliError, BirliError::DryRun},
    flags::FlagContext,
    io::IOContext,
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        hifitime::Epoch,
        mwalib::{CorrelatorContext, GeometricDelaysApplied},
        precession::{precess_time, PrecessionInfo},
        LatLngHeight, RADec,
    },
    passband_gains::{PFB_COTTER_2014_10KHZ, PFB_JAKE_2022_200HZ},
    Complex, PreprocessContext, VisSelection,
};
use cfg_if::cfg_if;
use clap::{arg, command, PossibleValue, ValueHint::FilePath};
use log::{debug, info, trace, warn};
use prettytable::{cell, format as prettyformat, row, table};
use std::{env, ffi::OsString, fmt::Debug};

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
pub fn display_build_info() {
    match GIT_HEAD_REF {
        Some(hr) => {
            let dirty = GIT_DIRTY.unwrap_or(false);
            info!(
                "Compiled on git commit hash: {}{}",
                GIT_COMMIT_HASH.unwrap(),
                if dirty { " (dirty)" } else { "" }
            );
            info!("            git head ref: {}", hr);
        }
        None => info!("Compiled on git commit hash: <no git info>"),
    }
    info!("            {}", BUILT_TIME_UTC);
    info!("         with compiler {}", RUSTC_VERSION);
    info!("");
}

/// Show info about Birli parameters
// TODO: fix too_many_arguments
#[allow(clippy::too_many_arguments)]
pub fn show_param_info(
    corr_ctx: &CorrelatorContext,
    prep_ctx: &PreprocessContext,
    flag_ctx: &FlagContext,
    vis_sel: &VisSelection,
    avg_time: usize,
    avg_freq: usize,
    num_timesteps_per_chunk: Option<usize>,
) {
    info!(
        "{} version {}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION"),
    );

    display_build_info();

    info!(
        "observation name:     {}",
        corr_ctx.metafits_context.obs_name
    );

    info!("Array position:       {}", &prep_ctx.array_pos);
    info!("Phase centre:         {}", &prep_ctx.phase_centre);
    let pointing_centre = RADec::from_mwalib_tile_pointing(&corr_ctx.metafits_context);
    if pointing_centre != prep_ctx.phase_centre {
        info!("Pointing centre:      {}", &pointing_centre);
    }

    let coarse_chan_flag_idxs: Vec<usize> = flag_ctx
        .coarse_chan_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    // TODO: actually display this.
    let _fine_chan_flag_idxs: Vec<usize> = flag_ctx
        .fine_chan_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    let timestep_flag_idxs: Vec<usize> = flag_ctx
        .timestep_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
    #[allow(clippy::needless_collect)]
    let baseline_flag_idxs: Vec<usize> = flag_ctx
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

    let (sched_start_date, sched_start_time, sched_start_mjd_s, sched_start_prec) = time_details(
        corr_ctx.metafits_context.sched_start_gps_time_ms,
        prep_ctx.phase_centre,
        prep_ctx.array_pos,
    );
    info!(
        "Scheduled start:      {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
        sched_start_date, sched_start_time,
        corr_ctx.metafits_context.sched_start_unix_time_ms as f64 / 1e3,
        corr_ctx.metafits_context.sched_start_gps_time_ms as f64 / 1e3,
        sched_start_mjd_s,
        sched_start_prec.lmst.to_degrees(),
        sched_start_prec.lmst_j2000.to_degrees(),
        sched_start_prec.array_latitude_j2000.to_degrees(),
    );
    let (sched_end_date, sched_end_time, sched_end_mjd_s, sched_end_prec) = time_details(
        corr_ctx.metafits_context.sched_end_gps_time_ms,
        prep_ctx.phase_centre,
        prep_ctx.array_pos,
    );
    info!(
        "Scheduled end:        {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
        sched_end_date, sched_end_time,
        corr_ctx.metafits_context.sched_end_unix_time_ms as f64 / 1e3,
        corr_ctx.metafits_context.sched_end_gps_time_ms as f64 / 1e3,
        sched_end_mjd_s,
        sched_end_prec.lmst.to_degrees(),
        sched_end_prec.lmst_j2000.to_degrees(),
        sched_end_prec.array_latitude_j2000.to_degrees(),
    );
    let int_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1e3;
    let sched_duration_s = corr_ctx.metafits_context.sched_duration_ms as f64 / 1e3;
    info!(
        "Scheduled duration:   {:.3}s = {:3} * {:.3}s",
        sched_duration_s,
        (sched_duration_s / int_time_s).ceil(),
        int_time_s
    );
    let quack_duration_s = corr_ctx.metafits_context.quack_time_duration_ms as f64 / 1e3;
    info!(
        "Quack duration:       {:.3}s = {:3} * {:.3}s",
        quack_duration_s,
        (quack_duration_s / int_time_s).ceil(),
        int_time_s
    );
    let num_avg_timesteps = (vis_sel.timestep_range.len() as f64 / avg_time as f64).ceil() as usize;
    let avg_int_time_s = int_time_s * avg_time as f64;
    info!(
        "Output duration:      {:.3}s = {:3} * {:.3}s{}",
        num_avg_timesteps as f64 * avg_int_time_s,
        num_avg_timesteps,
        avg_int_time_s,
        if avg_time != 1 {
            format!(" ({}x)", avg_time)
        } else {
            "".into()
        }
    );

    let total_bandwidth_mhz = corr_ctx.metafits_context.obs_bandwidth_hz as f64 / 1e6;
    let fine_chan_width_khz = corr_ctx.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

    info!(
        "Scheduled Bandwidth:  {:.3}MHz = {:3} * {:3} * {:.3}kHz",
        total_bandwidth_mhz,
        corr_ctx.metafits_context.num_metafits_coarse_chans,
        fine_chans_per_coarse,
        fine_chan_width_khz
    );

    let out_bandwidth_mhz =
        vis_sel.coarse_chan_range.len() as f64 * fine_chans_per_coarse as f64 * fine_chan_width_khz
            / 1e3;
    let num_avg_chans = (vis_sel.coarse_chan_range.len() as f64 * fine_chans_per_coarse as f64
        / avg_freq as f64)
        .ceil() as usize;
    let avg_fine_chan_width_khz = fine_chan_width_khz * avg_freq as f64;
    info!(
        "Output Bandwidth:     {:.3}MHz = {:9} * {:.3}kHz{}",
        out_bandwidth_mhz,
        num_avg_chans,
        avg_fine_chan_width_khz,
        if avg_freq != 1 {
            format!(" ({}x)", avg_freq)
        } else {
            "".into()
        }
    );

    let first_epoch = Epoch::from_gpst_seconds(corr_ctx.timesteps[0].gps_time_ms as f64 / 1e3);
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

    let provided_timestep_indices = corr_ctx.provided_timestep_indices.clone();
    let common_timestep_indices = corr_ctx.common_timestep_indices.clone();
    let common_good_timestep_indices = corr_ctx.common_good_timestep_indices.clone();
    for (timestep_idx, timestep) in corr_ctx.timesteps.iter().enumerate() {
        let provided = provided_timestep_indices.contains(&timestep_idx);
        let selected = vis_sel.timestep_range.contains(&timestep_idx);
        let common = common_timestep_indices.contains(&timestep_idx);
        let good = common_good_timestep_indices.contains(&timestep_idx);
        let flagged = timestep_flag_idxs.contains(&timestep_idx);

        let (_, time, ..) = time_details(
            timestep.gps_time_ms,
            prep_ctx.phase_centre,
            prep_ctx.array_pos,
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

    info!(
        "Timestep details (all={}, provided={}, common={}, good={}, select={}, flag={}):{}",
        corr_ctx.num_timesteps,
        corr_ctx.num_provided_timesteps,
        corr_ctx.num_common_timesteps,
        corr_ctx.num_common_good_timesteps,
        vis_sel.timestep_range.len(),
        timestep_flag_idxs.len(),
        if show_timestep_table {
            format!("\n{}", timestep_table)
        } else {
            "".into()
        }
    );
    if !show_timestep_table {
        info!("-> provided:    {:?}", corr_ctx.provided_timestep_indices);
        info!("-> common:      {:?}", corr_ctx.common_timestep_indices);
        info!(
            "-> common good: {:?}",
            corr_ctx.common_good_timestep_indices
        );
        info!("-> selected:    {:?}", vis_sel.timestep_range);
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
    let provided_coarse_chan_indices = corr_ctx.provided_coarse_chan_indices.clone();
    let common_coarse_chan_indices = corr_ctx.common_coarse_chan_indices.clone();
    let common_good_coarse_chan_indices = corr_ctx.common_good_coarse_chan_indices.clone();
    for (chan_idx, chan) in corr_ctx.coarse_chans.iter().enumerate() {
        let provided = provided_coarse_chan_indices.contains(&chan_idx);
        let selected = vis_sel.coarse_chan_range.contains(&chan_idx);
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

    info!(
        "Coarse channel details (metafits={}, provided={}, common={}, good={}, select={}, flag={}):{}",
        corr_ctx.num_coarse_chans,
        corr_ctx.num_provided_coarse_chans,
        corr_ctx.num_common_coarse_chans,
        corr_ctx.num_common_good_coarse_chans,
        vis_sel.coarse_chan_range.len(),
        coarse_chan_flag_idxs.len(),
        if show_coarse_chan_table { format!("\n{}", coarse_chan_table) } else { "".into() }
    );

    if !show_coarse_chan_table {
        info!(
            "-> provided:    {:?}",
            corr_ctx.provided_coarse_chan_indices
        );
        info!("-> common:      {:?}", corr_ctx.common_coarse_chan_indices);
        info!(
            "-> common good: {:?}",
            corr_ctx.common_good_coarse_chan_indices
        );
        info!("-> selected:    {:?}", vis_sel.coarse_chan_range);
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

    for (ant_idx, ant) in corr_ctx.metafits_context.antennas.iter().enumerate() {
        let flagged = *flag_ctx.antenna_flags.get(ant_idx).unwrap_or(&false);
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
        corr_ctx.metafits_context.num_ants,
        flag_ctx
            .antenna_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .count(),
        format!("\n{}", ant_table)
    );

    // let show_baseline_table = false;

    info!(
        "Baseline Details (all={}, auto={}, select={}, flag={}):",
        corr_ctx.metafits_context.num_baselines,
        corr_ctx.metafits_context.num_ants,
        vis_sel.baseline_idxs.len(),
        baseline_flag_idxs.len(),
    );

    // if !show_baseline_table {
    //     info!("-> selected:    {:?}", vis_sel.baseline_idxs);
    //     info!("-> flags:    {:?}", baseline_flag_idxs);
    // }

    // TODO: show free memory with https://docs.rs/sys-info/latest/sys_info/fn.mem_info.html

    let num_sel_timesteps = vis_sel.timestep_range.len();
    let num_sel_chans = vis_sel.coarse_chan_range.len() * fine_chans_per_coarse;
    let num_sel_baselines = vis_sel.baseline_idxs.len();
    let num_sel_pols = corr_ctx.metafits_context.num_visibility_pols;
    let mem_per_timestep_gib = (num_sel_chans
        * num_sel_baselines
        * num_sel_pols
        * (std::mem::size_of::<Complex<f32>>()
            + std::mem::size_of::<f32>()
            + std::mem::size_of::<bool>())) as f64
        / 1024.0_f64.powi(3);

    info!(
        "Estimated memory usage per timestep =           {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
        num_sel_chans,
        num_sel_baselines,
        num_sel_pols,
        std::mem::size_of::<Complex<f32>>(),
        std::mem::size_of::<f32>(),
        std::mem::size_of::<bool>(),
        mem_per_timestep_gib,
    );

    if let Some(num_timesteps) = num_timesteps_per_chunk {
        info!("Estimated memory per chunk          = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
            num_timesteps,
            num_sel_chans,
            num_sel_baselines,
            num_sel_pols,
            std::mem::size_of::<Complex<f32>>(),
            std::mem::size_of::<f32>(),
            std::mem::size_of::<bool>(),
            mem_per_timestep_gib * num_timesteps as f64,
        );
    }

    info!("Estimated memory selected           = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
        num_sel_timesteps,
        num_sel_chans,
        num_sel_baselines,
        num_sel_pols,
        std::mem::size_of::<Complex<f32>>(),
        std::mem::size_of::<f32>(),
        std::mem::size_of::<bool>(),
        mem_per_timestep_gib * num_sel_timesteps as f64,
    );

    let avg_mem_per_timestep_gib = (num_avg_chans
        * num_sel_baselines
        * num_sel_pols
        * (std::mem::size_of::<Complex<f32>>()
            + std::mem::size_of::<f32>()
            + std::mem::size_of::<bool>())) as f64
        / 1024.0_f64.powi(3);

    info!("Estimated output size               = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>) = {:7.02} GiB",
        num_avg_timesteps,
        num_avg_chans,
        num_sel_baselines,
        num_sel_pols,
        std::mem::size_of::<Complex<f32>>(),
        std::mem::size_of::<f32>(),
        std::mem::size_of::<bool>(),
        avg_mem_per_timestep_gib * num_avg_timesteps as f64,
    );
}

/// Parse an iterator of arguments into a BirliContext.
///
/// # Errors
///
/// mostly just parsing stuff in here.
pub fn parse_args<I, T>(args: I) -> Result<BirliContext, BirliError>
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
            arg!(--"sel-time" "[WIP] Timestep index range (inclusive) to select")
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
            // TODO: rename to antennas
            // -> antennae
            arg!(--"no-flag-metafits" "[WIP] Ignore antenna flags in metafits")
                .help_heading("FLAGGING"),
            arg!(--"flag-antennae" <ANTS>... "[WIP] Flag antenna indices")
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

    // base command line matches
    let matches = app.get_matches_from(args);
    trace!("arg matches:\n{:?}", &matches);

    // // //
    // IO //
    // // //

    // TODO: proper error handling instead of .expect()
    let io_ctx = IOContext {
        metafits_in: matches
            .value_of("metafits")
            .expect("--metafits must be a valid path")
            .into(),
        gpufits_in: matches
            .values_of("fits_paths")
            .expect("--")
            .map(|s| s.into())
            .collect(),
        aocalsols_in: matches.value_of("apply-di-cal").map(|s| s.into()),
        uvfits_out: matches.value_of("uvfits-out").map(|s| s.into()),
        ms_out: matches.value_of("ms-out").map(|s| s.into()),
        flag_template: matches.value_of("flag-template").map(|s| s.into()),
    };

    let corr_ctx = io_ctx.get_corr_ctx().expect("unable to get mwalib context");
    debug!("mwalib correlator context:\n{}", &corr_ctx);

    // ////////// //
    // Selections //
    // ////////// //

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    if let Some(mut values) = matches.values_of("sel-time") {
        // TODO: custom error types
        if let (Some(from), Some(to)) = (values.next(), values.next()) {
            let from = from.parse::<usize>().expect("cannot parse --sel-time from");
            assert!(
                from < corr_ctx.num_timesteps,
                "invalid --sel-time from {}. must be < num_timesteps ({})",
                from,
                corr_ctx.num_timesteps
            );
            let to = to.parse::<usize>().expect("cannot parse --sel-time to");
            assert!(
                to >= from && to < corr_ctx.num_timesteps,
                "invalid --sel-time from {} to {}, must be < num_timesteps ({})",
                from,
                to,
                corr_ctx.num_timesteps
            );
            vis_sel.timestep_range = from..to + 1
        } else {
            panic!("invalid --sel-time <from> <to>, two values must be provided");
        }
    };

    // //////// //
    // Flagging //
    // //////// //

    let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);

    // Timesteps
    if let Some(values) = matches.values_of("flag-times") {
        values
            .map(|value| {
                // TODO: custom error types
                if let Ok(timestep_idx) = value.parse::<usize>() {
                    flag_ctx.timestep_flags[timestep_idx] = true;
                } else {
                    panic!("unable to parse timestep value: {}", value);
                }
            })
            .collect()
    };

    // TODO: init and end steps
    // let mut init_steps: usize = 0;
    // let mut end_steps: usize = 0;
    // if let Some(count_str) = matches.value_of("flag-init-steps") {
    //     init_steps = count_str.parse::<usize>().unwrap();
    //     info!("Flagging {} initial timesteps", init_steps);
    // }
    // if let Some(seconds_str) = matches.value_of("flag-init") {
    //     let init_seconds = seconds_str.parse::<f64>().unwrap();
    //     // init_steps = todo!();
    // }

    // coarse channels
    if let Some(coarse_chans) = matches.values_of("flag-coarse-chans") {
        for value in coarse_chans {
            // TODO: custom error types
            if let Ok(coarse_chan_idx) = value.parse::<usize>() {
                flag_ctx.coarse_chan_flags[coarse_chan_idx] = true;
            } else {
                panic!("unable to parse coarse chan value: {}", value);
            }
        }
    }

    // fine channels
    if let Some(fine_chans) = matches.values_of("flag-fine-chans") {
        for value in fine_chans {
            // TODO: custom error types
            if let Ok(fine_chan_idx) = value.parse::<usize>() {
                flag_ctx.fine_chan_flags[fine_chan_idx] = true;
            } else {
                panic!("unable to parse fine_chan value: {}", value);
            }
        }
    }

    // Antennas
    let ignore_metafits = matches.is_present("no-flag-metafits");
    if ignore_metafits {
        info!("Ignoring antenna flags from metafits.");
        // set antenna flags to all false
        flag_ctx.antenna_flags = vec![false; flag_ctx.antenna_flags.len()];
    }

    if let Some(antennae) = matches.values_of("flag-antennae") {
        for value in antennae {
            // TODO: custom error types
            if let Ok(antenna_idx) = value.parse::<usize>() {
                flag_ctx.antenna_flags[antenna_idx] = true;
            } else {
                panic!("unable to parse antenna value: {}", value);
            }
        }
    }

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
        "flag-antennae",
        "sel-time",
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

    show_param_info(
        &corr_ctx,
        &prep_ctx,
        &flag_ctx,
        &vis_sel,
        avg_time,
        avg_freq,
        num_timesteps_per_chunk,
    );

    prep_ctx.log_info();

    if matches.is_present("dry-run") {
        return Err(DryRun {});
    }

    Ok(BirliContext {
        corr_ctx,
        prep_ctx,
        vis_sel,
        flag_ctx,
        io_ctx,
        avg_time,
        avg_freq,
        num_timesteps_per_chunk,
    })
}
