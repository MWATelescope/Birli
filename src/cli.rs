//! Command line interface.

use crate::{
    calibration::apply_di_calsol,
    context_to_jones_array, correct_cable_lengths, correct_geometry,
    corrections::correct_digital_gains,
    flags::{
        add_dimension, flag_to_weight_array, get_baseline_flags, get_coarse_chan_flags,
        get_coarse_chan_range, get_timestep_flags, get_timestep_range, get_weight_factor,
    },
    get_antenna_flags, init_flag_array,
    io::{aocal::AOCalSols, WriteableVis},
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib::{CorrelatorContext, GeometricDelaysApplied},
        LatLngHeight,
    },
    write_flags, Axis, Complex, UvfitsWriter,
};
use cfg_if::cfg_if;
use clap::{app_from_crate, arg, AppSettings, ValueHint::FilePath};
use itertools::Itertools;
use log::{debug, info, trace, warn};
use marlu::{
    hifitime::Epoch,
    io::{ms::MeasurementSetWriter, VisWritable},
    precession::{precess_time, PrecessionInfo},
    RADec,
};
use prettytable::{cell, format as prettyformat, row, table};
use std::{collections::HashMap, env, ffi::OsString, fmt::Debug, ops::Range, time::Duration};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use crate::{
            flags::flag_jones_array_existing,
        };
        use aoflagger_sys::{cxx_aoflagger_new};
    }
}

// Add build-time information from the "built" crate.
include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// stolen from hyperdrive
/// Write many info-level log lines of how this executable was compiled.
fn display_build_info() {
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

// TODO: fix too_many_arguments
#[allow(clippy::too_many_arguments)]
fn show_param_info(
    context: &CorrelatorContext,
    array_pos: LatLngHeight,
    phase_centre: RADec,
    coarse_chan_range: &Range<usize>,
    timestep_range: &Range<usize>,
    baseline_idxs: &[usize],
    coarse_chan_flags: &[bool],
    fine_chan_flags: &[bool],
    timestep_flags: &[bool],
    antenna_flags: &[bool],
    baseline_flags: &[bool],
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
        context.metafits_context.obs_name
    );

    info!("Array position:       {}", &array_pos);
    info!("Phase centre:         {}", &phase_centre);
    let pointing_centre = RADec::from_mwalib_tile_pointing(&context.metafits_context);
    if pointing_centre != phase_centre {
        info!("Pointing centre:      {}", &pointing_centre);
    }

    let coarse_chan_flag_idxs: Vec<usize> = coarse_chan_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    // TODO: actually display this.
    let _fine_chan_flag_idxs: Vec<usize> = fine_chan_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    let timestep_flag_idxs: Vec<usize> = timestep_flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
        .collect();
    #[allow(clippy::needless_collect)]
    let baseline_flag_idxs: Vec<usize> = baseline_flags
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
        context.metafits_context.sched_start_gps_time_ms,
        phase_centre,
        array_pos,
    );
    info!(
        "Scheduled start:      {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
        sched_start_date, sched_start_time,
        context.metafits_context.sched_start_unix_time_ms as f64 / 1e3,
        context.metafits_context.sched_start_gps_time_ms as f64 / 1e3,
        sched_start_mjd_s,
        sched_start_prec.lmst.to_degrees(),
        sched_start_prec.lmst_j2000.to_degrees(),
        sched_start_prec.array_latitude_j2000.to_degrees(),
    );
    let (sched_end_date, sched_end_time, sched_end_mjd_s, sched_end_prec) = time_details(
        context.metafits_context.sched_end_gps_time_ms,
        phase_centre,
        array_pos,
    );
    info!(
        "Scheduled end:        {} {} UTC, unix={:.3}, gps={:.3}, mjd={:.3}, lmst={:7.4}°, lmst2k={:7.4}°, lat2k={:7.4}°",
        sched_end_date, sched_end_time,
        context.metafits_context.sched_end_unix_time_ms as f64 / 1e3,
        context.metafits_context.sched_end_gps_time_ms as f64 / 1e3,
        sched_end_mjd_s,
        sched_end_prec.lmst.to_degrees(),
        sched_end_prec.lmst_j2000.to_degrees(),
        sched_end_prec.array_latitude_j2000.to_degrees(),
    );
    let int_time_s = context.metafits_context.corr_int_time_ms as f64 / 1e3;
    let sched_duration_s = context.metafits_context.sched_duration_ms as f64 / 1e3;
    info!(
        "Scheduled duration:   {:.3}s = {:3} * {:.3}s",
        sched_duration_s,
        (sched_duration_s / int_time_s).ceil(),
        int_time_s
    );
    let quack_duration_s = context.metafits_context.quack_time_duration_ms as f64 / 1e3;
    info!(
        "Quack duration:       {:.3}s = {:3} * {:.3}s",
        quack_duration_s,
        (quack_duration_s / int_time_s).ceil(),
        int_time_s
    );
    let num_avg_timesteps = (timestep_range.len() as f64 / avg_time as f64).ceil() as usize;
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

    let total_bandwidth_mhz = context.metafits_context.obs_bandwidth_hz as f64 / 1e6;
    let fine_chan_width_khz = context.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

    info!(
        "Scheduled Bandwidth:  {:.3}MHz = {:3} * {:3} * {:.3}kHz",
        total_bandwidth_mhz,
        context.metafits_context.num_metafits_coarse_chans,
        fine_chans_per_coarse,
        fine_chan_width_khz
    );

    let out_bandwidth_mhz =
        coarse_chan_range.len() as f64 * fine_chans_per_coarse as f64 * fine_chan_width_khz / 1e3;
    let num_avg_chans = (coarse_chan_range.len() as f64 * fine_chans_per_coarse as f64
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

    let first_epoch = Epoch::from_gpst_seconds(context.timesteps[0].gps_time_ms as f64 / 1e3);
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

    let provided_timestep_indices = context.provided_timestep_indices.clone();
    let common_timestep_indices = context.common_timestep_indices.clone();
    let common_good_timestep_indices = context.common_good_timestep_indices.clone();
    for (timestep_idx, timestep) in context.timesteps.iter().enumerate() {
        let provided = provided_timestep_indices.contains(&timestep_idx);
        let selected = timestep_range.contains(&timestep_idx);
        let common = common_timestep_indices.contains(&timestep_idx);
        let good = common_good_timestep_indices.contains(&timestep_idx);
        let flagged = timestep_flag_idxs.contains(&timestep_idx);

        let (_, time, ..) = time_details(timestep.gps_time_ms, phase_centre, array_pos);
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
        context.num_timesteps,
        context.num_provided_timesteps,
        context.num_common_timesteps,
        context.num_common_good_timesteps,
        timestep_range.len(),
        timestep_flag_idxs.len(),
        if show_timestep_table {
            format!("\n{}", timestep_table)
        } else {
            "".into()
        }
    );
    if !show_timestep_table {
        info!("-> provided:    {:?}", context.provided_timestep_indices);
        info!("-> common:      {:?}", context.common_timestep_indices);
        info!("-> common good: {:?}", context.common_good_timestep_indices);
        info!("-> selected:    {:?}", timestep_range);
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
    let provided_coarse_chan_indices = context.provided_coarse_chan_indices.clone();
    let common_coarse_chan_indices = context.common_coarse_chan_indices.clone();
    let common_good_coarse_chan_indices = context.common_good_coarse_chan_indices.clone();
    for (chan_idx, chan) in context.coarse_chans.iter().enumerate() {
        let provided = provided_coarse_chan_indices.contains(&chan_idx);
        let selected = coarse_chan_range.contains(&chan_idx);
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
        context.num_coarse_chans,
        context.num_provided_coarse_chans,
        context.num_common_coarse_chans,
        context.num_common_good_coarse_chans,
        coarse_chan_range.len(),
        coarse_chan_flag_idxs.len(),
        if show_coarse_chan_table { format!("\n{}", coarse_chan_table) } else { "".into() }
    );

    if !show_coarse_chan_table {
        info!("-> provided:    {:?}", context.provided_coarse_chan_indices);
        info!("-> common:      {:?}", context.common_coarse_chan_indices);
        info!(
            "-> common good: {:?}",
            context.common_good_coarse_chan_indices
        );
        info!("-> selected:    {:?}", coarse_chan_range);
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

    for (ant_idx, ant) in context.metafits_context.antennas.iter().enumerate() {
        let flagged = *antenna_flags.get(ant_idx).unwrap_or(&false);
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
        context.metafits_context.num_ants,
        antenna_flags
            .iter()
            .enumerate()
            .filter_map(|(idx, &flag)| if flag { Some(idx) } else { None })
            .count(),
        format!("\n{}", ant_table)
    );

    // let show_baseline_table = false;

    info!(
        "Baseline Details (all={}, auto={}, select={}, flag={}):",
        context.metafits_context.num_baselines,
        context.metafits_context.num_ants,
        baseline_idxs.len(),
        baseline_flag_idxs.len(),
    );

    // if !show_baseline_table {
    //     info!("-> selected:    {:?}", baseline_idxs);
    //     info!("-> flags:    {:?}", baseline_flag_idxs);
    // }

    // TODO: show:
    // - estimated memory consumption
    // - free memory with https://docs.rs/sys-info/latest/sys_info/fn.mem_info.html

    let num_sel_timesteps = timestep_range.len();
    let num_sel_chans = coarse_chan_range.len() * fine_chans_per_coarse;
    let num_sel_baselines = baseline_idxs.len();
    let num_sel_pols = context.metafits_context.num_visibility_pols;
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

macro_rules! with_increment_duration {
    ($durs:expr, $name:literal, $($s:stmt);+ $(;)?) => {
        {
            let _now = std::time::Instant::now();
            let _res = {
                $(
                    $s
                )*
            };
            *$durs.entry($name).or_insert(Duration::default()) += _now.elapsed();
            _res
        }
    };
}

/// Main Birli preprocessing function
pub fn main_with_args<I, T>(args: I)
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
    I: Debug,
{
    debug!("args:\n{:?}", &args);

    // TODO: fix this
    #[allow(unused_mut)]
    let mut app = app_from_crate!()
        .setting(AppSettings::SubcommandPrecedenceOverArg)
        .setting(AppSettings::ArgRequiredElseHelp)
        .unset_setting(AppSettings::NextLineHelp)
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
            // TODO make this work the same way as rust ranges. start <= x < end
            arg!(--"sel-time" "[WIP] Timestep index range (inclusive) to select")
                .help_heading("SELECTION")
                .value_names(&["MIN", "MAX"])
                .required(false),
            arg!(--"no-sel-flagged-ants" "[WIP] Deselect flagged antennas")
                .help_heading("SELECTION")
                .required(false),
            arg!(--"no-sel-autos" "[WIP] Deselect autocorrelations")
                .help_heading("SELECTION")
                .required(false),

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
                .help_heading("FLAGGING")
                .required(false),
            arg!(--"no-flag-dc" "[WIP] Do not flag DC centre chans")
                .help_heading("FLAGGING")
                .required(false),
            // -> antennae
            arg!(--"no-flag-metafits" "[WIP] Ignore antenna flags in metafits")
                .help_heading("FLAGGING")
                .required(false),
            arg!(--"flag-antennae" <ANTS>... "[WIP] Flag antenna indices")
                .help_heading("FLAGGING")
                .multiple_values(true)
                .required(false),
            // -> baselines
            arg!(--"flag-autos" "[WIP] Flag auto correlations")
                .help_heading("FLAGGING")
                .required(false),

            // corrections
            arg!(--"no-cable-delay" "Do not perform cable length corrections")
                .help_heading("CORRECTION"),
            arg!(--"no-geometric-delay" "Do not perform geometric corrections")
                .help_heading("CORRECTION")
                .alias("no-geom"),
            arg!(--"no-digital-gains" "Do not perform digital gains corrections")
                .help_heading("CORRECTION"),
            arg!(--"passband" <PATH> "[WIP] Apply passband corrections from <PATH>")
                .required(false)
                .value_hint(FilePath)
                .help_heading("CORRECTION"),

            // calibration
            arg!(--"apply-di-cal" <PATH> "Apply DI calibration solutions to the data before averaging")
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
        } else {
            // add the --no-rfi flag anyway, but hidden.
            app = app.args(&[
                arg!(--"no-rfi" "Does nothing without the aoflagger feature")
            ]);
        }
    };

    // base command line matches
    let matches = app.get_matches_from(args);
    trace!("arg matches:\n{:?}", &matches);

    let metafits_path = matches
        .value_of("metafits")
        .expect("--metafits must be a valid path");
    let fits_paths: Vec<&str> = matches.values_of("fits_paths").expect("--").collect();

    let context =
        CorrelatorContext::new(&metafits_path, &fits_paths).expect("unable to get mwalib context");
    debug!("mwalib correlator context:\n{}", &context);
    let sel_coarse_chan_range = get_coarse_chan_range(&context).unwrap();
    let mut coarse_chan_flags = get_coarse_chan_flags(&context);
    let sel_timestep_range = match matches.values_of("sel-time") {
        Some(mut values) => {
            if let (Some(from), Some(to)) = (values.next(), values.next()) {
                let from = from.parse::<usize>().expect("cannot parse --sel-time from");
                assert!(
                    from < context.num_timesteps,
                    "invalid --sel-time from {}. must be < num_timesteps ({})",
                    from,
                    context.num_timesteps
                );
                let to = to.parse::<usize>().expect("cannot parse --sel-time to");
                assert!(
                    to >= from && to < context.num_timesteps,
                    "invalid --sel-time from {} to {}, must be < num_timesteps ({})",
                    from,
                    to,
                    context.num_timesteps
                );
                from..to + 1
            } else {
                panic!("invalid --sel-time <from> <to>, two values must be provided");
            }
        }
        _ => get_timestep_range(&context).unwrap(),
    };

    let mut timestep_flags = get_timestep_flags(&context);
    let mut antenna_flags = get_antenna_flags(&context);
    let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

    let mut fine_chan_flags = vec![false; fine_chans_per_coarse];

    let array_pos = if matches.is_present("emulate-cotter") {
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

    let phase_centre = match (
        matches.values_of("phase-centre"),
        matches.is_present("pointing-centre"),
    ) {
        (Some(mut values), _) => {
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
        (_, true) => RADec::from_mwalib_tile_pointing(&context.metafits_context),
        _ => RADec::from_mwalib_phase_or_pointing(&context.metafits_context),
    };

    // /////////////// //
    // Manual flagging //
    // /////////////// //

    // coarse channels
    if let Some(coarse_chans) = matches.values_of("flag-coarse-chans") {
        for value in coarse_chans {
            if let Ok(coarse_chan_idx) = value.parse::<usize>() {
                coarse_chan_flags[coarse_chan_idx] = true;
            } else {
                panic!("unable to parse coarse chan value: {}", value);
            }
        }
    }

    // fine channels
    if let Some(fine_chans) = matches.values_of("flag-fine-chans") {
        for value in fine_chans {
            if let Ok(fine_chan_idx) = value.parse::<usize>() {
                fine_chan_flags[fine_chan_idx] = true;
            } else {
                panic!("unable to parse fine_chan value: {}", value);
            }
        }
    }

    // time
    // TODO
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
    if let Some(timesteps) = matches.values_of("flag-times") {
        for value in timesteps {
            if let Ok(timestep_idx) = value.parse::<usize>() {
                timestep_flags[timestep_idx] = true;
            } else {
                panic!("unable to parse timestep value: {}", value);
            }
        }
    }

    // antennae
    // TODO
    let ignore_metafits = matches.is_present("no-flag-metafits");
    if ignore_metafits {
        info!("Ignoring antenna flags from metafits.");
        // set antenna flags to all false
        antenna_flags = vec![false; antenna_flags.len()];
    }

    if let Some(antennae) = matches.values_of("flag-antennae") {
        for value in antennae {
            if let Ok(antenna_idx) = value.parse::<usize>() {
                antenna_flags[antenna_idx] = true;
            } else {
                panic!("unable to parse antenna value: {}", value);
            }
        }
    }
    // } else {
    //     let init_seconds = context.metafits_context.quack_time_duration_ms as f64 / 1e3;
    //     let edge_width_khz = 40.0;
    //     info!("Using default flagging parameters. {} seconds, {} kHz edges", init_seconds, edge_width_khz);
    //     // init_steps = todo!();
    // }

    // ///////// //
    // Averaging //
    // ///////// //

    let int_time_s = context.metafits_context.corr_int_time_ms as f64 / 1e3;

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

    let fine_chan_width_khz = context.metafits_context.corr_fine_chan_width_hz as f64 / 1e3;

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

    let num_sel_timesteps = sel_timestep_range.len();
    let num_sel_chans = sel_coarse_chan_range.len() * fine_chans_per_coarse;
    let num_sel_baselines = baseline_idxs.len();
    let num_sel_pols = context.metafits_context.num_visibility_pols;
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
            panic!("you can't use --time-chunk and --max-memory at the same time");
        }
        (Some(steps_str), None) => {
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

    let baseline_flags = get_baseline_flags(&context, &antenna_flags);

    // ///////// //
    // Show info //
    // ///////// //

    show_param_info(
        &context,
        array_pos,
        phase_centre,
        &sel_coarse_chan_range,
        &sel_timestep_range,
        &baseline_idxs,
        &coarse_chan_flags,
        &fine_chan_flags,
        &timestep_flags,
        &antenna_flags,
        &baseline_flags,
        avg_time,
        avg_freq,
        num_timesteps_per_chunk,
    );

    if matches.is_present("dry-run") {
        info!("Dry run. No files will be written.");
        return;
    }

    let draw_progress = !matches.is_present("no-draw-progress");

    for unimplemented_option in &[
        "flag-init",
        "flag-init-steps",
        "flag-autos",
        "flag-end",
        "flag-end-steps",
        "flag-edge-width",
        "flag-edge-chans",
        "flag-fine-chans",
        "flag-dc",
        "no-flag-dc",
        "passband",
        "no-sel-autos",
        "no-sel-flagged-ants",
    ] {
        if matches.is_present(unimplemented_option) {
            panic!("option not yet implemented: --{}", unimplemented_option);
        }
    }

    for untested_option in &[
        "flag-times",
        "flag-coarse-chans",
        "flag-fine-chans",
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

    let sel_baselines = baseline_idxs
        .iter()
        .map(|&idx| {
            let baseline = &context.metafits_context.baselines[idx];
            (baseline.ant1_index, baseline.ant2_index)
        })
        .collect::<Vec<_>>();

    // used to time large operations
    let mut durations = HashMap::<&str, Duration>::new();

    let mut uvfits_writer = matches.value_of("uvfits-out").map(|uvfits_out| {
        with_increment_duration!(durations, "init", {
            UvfitsWriter::from_mwalib(
                uvfits_out,
                &context,
                &sel_timestep_range,
                &sel_coarse_chan_range,
                &baseline_idxs,
                Some(array_pos),
                Some(phase_centre),
                avg_time,
                avg_freq,
            )
            .expect("couldn't initialise uvfits writer")
        })
    });
    let mut ms_writer = matches.value_of("ms-out").map(|ms_out| {
        let writer = MeasurementSetWriter::new(ms_out, phase_centre, Some(array_pos));
        with_increment_duration!(durations, "init", {
            writer
                .initialize_from_mwalib(
                    &context,
                    &sel_timestep_range,
                    &sel_coarse_chan_range,
                    &baseline_idxs,
                    avg_time,
                    avg_freq,
                )
                .unwrap();
        });
        writer
    });

    let chunk_size = if let Some(steps) = num_timesteps_per_chunk {
        steps
    } else {
        sel_timestep_range.len()
    };
    for mut timestep_chunk in &sel_timestep_range.clone().chunks(chunk_size) {
        let chunk_first_timestep = timestep_chunk.next().unwrap();
        let chunk_last_timestep = timestep_chunk.last().unwrap_or(chunk_first_timestep);
        let chunk_timestep_range: Range<usize> = chunk_first_timestep..chunk_last_timestep + 1;
        if num_timesteps_per_chunk.is_some() {
            info!(
                "processing timestep chunk {:?} of {:?} % {}",
                chunk_timestep_range, sel_timestep_range, chunk_size
            );
        }
        let flag_array = init_flag_array(
            &context,
            &chunk_timestep_range,
            &sel_coarse_chan_range,
            Some(&timestep_flags),
            Some(&coarse_chan_flags),
            Some(&fine_chan_flags),
            Some(&baseline_flags),
        );

        #[allow(unused_mut)]
        let (mut jones_array, mut flag_array) = with_increment_duration!(
            durations,
            "read",
            context_to_jones_array(
                &context,
                &chunk_timestep_range,
                &sel_coarse_chan_range,
                Some(flag_array),
                draw_progress,
            )
            .unwrap()
        );

        // perform cable delays if user has not disabled it, and they haven't aleady beeen applied.

        let no_cable_delays = matches.is_present("no-cable-delay");
        let cable_delays_applied = context.metafits_context.cable_delays_applied;
        if !cable_delays_applied && !no_cable_delays {
            info!(
                "Applying cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
            with_increment_duration!(
                durations,
                "correct",
                correct_cable_lengths(&context, &mut jones_array, &sel_coarse_chan_range, false)
            );
        } else {
            info!(
                "Skipping cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
        }

        // perform coarse channel gain and passband corrections.
        if !matches.is_present("no-digital-gains") {
            with_increment_duration!(
                durations,
                "correct",
                correct_digital_gains(
                    &context,
                    &mut jones_array,
                    &sel_coarse_chan_range,
                    &sel_baselines,
                )
                .expect("couldn't apply digital gains")
            );
        }

        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                if !matches.is_present("no-rfi") {
                    let aoflagger = unsafe { cxx_aoflagger_new() };
                    let default_strategy_filename = aoflagger.FindStrategyFileMWA();
                    let strategy_filename = matches.value_of("aoflagger-strategy").unwrap_or(&default_strategy_filename);
                    info!("flagging with strategy {}", strategy_filename);
                    flag_array = with_increment_duration!(durations,
                        "flag",
                        flag_jones_array_existing(
                            &aoflagger,
                            strategy_filename,
                            &jones_array,
                            Some(flag_array),
                            true,
                            draw_progress,
                        )
                    );
                } else {
                    info!("skipped aoflagger");
                }
            }
        }

        // perform geometric delays if user has not disabled it, and they haven't aleady beeen applied.
        let no_geometric_delays = matches.is_present("no-geometric-delay");
        let geometric_delays_applied = context.metafits_context.geometric_delays_applied;
        match (geometric_delays_applied, no_geometric_delays) {
            (GeometricDelaysApplied::No, false) => {
                info!(
                    "Applying geometric delays. applied: {:?}, desired: {}",
                    geometric_delays_applied, !no_geometric_delays
                );
                with_increment_duration!(
                    durations,
                    "correct",
                    correct_geometry(
                        &context,
                        &mut jones_array,
                        &chunk_timestep_range,
                        &sel_coarse_chan_range,
                        Some(array_pos),
                        Some(phase_centre),
                        false,
                    )
                );
            }
            (..) => {
                info!(
                    "Skipping geometric delays. applied: {:?}, desired: {}",
                    geometric_delays_applied, !no_geometric_delays
                );
            }
        };

        // output flags (before averaging)
        if let Some(flag_template) = matches.value_of("flag-template") {
            with_increment_duration!(
                durations,
                "write",
                write_flags(&context, &flag_array, flag_template, &sel_coarse_chan_range)
                    .expect("unable to write flags")
            );
        }

        // generate weights
        let num_pols = context.metafits_context.num_visibility_pols;
        let weight_factor = get_weight_factor(&context);
        let mut weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        // let marlu_context = MarluVisContext::from_mwalib(
        //     &context,
        //     &chunk_timestep_range,
        //     &coarse_chan_range,
        //     &baseline_idxs,
        //     avg_time,
        //     avg_freq,
        // );

        if let Some(calsol_file) = matches.value_of("apply-di-cal") {
            let calsols = with_increment_duration!(
                durations,
                "read",
                AOCalSols::read_andre_binary(calsol_file).unwrap()
            );
            if calsols.di_jones.dim().0 != 1 {
                panic!("only 1 timeblock supported for calsols");
            }

            with_increment_duration!(
                durations,
                "calibrate",
                apply_di_calsol(
                    calsols.di_jones.index_axis(Axis(0), 0),
                    jones_array.view_mut(),
                    weight_array.view_mut(),
                    flag_array.view_mut(),
                    &sel_baselines,
                )
                .unwrap()
            );
        }

        // TODO: nothing actually uses the pol axis for flags and weights, so rip it out.
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
                        &context,
                        &chunk_timestep_range,
                        &sel_coarse_chan_range,
                        &baseline_idxs,
                        avg_time,
                        avg_freq,
                        draw_progress,
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
                        &context,
                        &chunk_timestep_range,
                        &sel_coarse_chan_range,
                        &baseline_idxs,
                        avg_time,
                        avg_freq,
                        draw_progress,
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
                .write_ants_from_mwalib(&context.metafits_context)
                .expect("couldn't write antenna table to uvfits")
        );
    }

    for (name, duration) in durations {
        info!("{} duration: {:?}", name, duration);
    }
}

#[cfg(test)]
mod tests {
    use std::env;

    use super::main_with_args;

    #[test]
    #[ignore = "flaky"]
    fn main_with_version_doesnt_crash() {
        main_with_args(&["birli", "--version"]);
    }

    #[test]
    #[ignore = "flaky"]
    fn forked_main_with_version_prints_version() {
        let pkg_name = env!("CARGO_PKG_NAME");
        let pkg_version = env!("CARGO_PKG_VERSION");
        assert_cli::Assert::main_binary()
            .with_args(&["--version"])
            .succeeds()
            .stdout()
            .contains(format!("{} {}\n", pkg_name, pkg_version).as_str())
            .unwrap();
    }
}
