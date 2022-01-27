use birli::{
    flags::{
        expand_flag_array, flag_to_weight_array, get_baseline_flags, get_coarse_chan_flags,
        get_coarse_chan_range, get_timestep_flags, get_timestep_range, get_weight_factor,
    },
    io::WriteableVis,
    UvfitsWriter,
};
use cfg_if::cfg_if;
use clap::{app_from_crate, arg, AppSettings, ValueHint::FilePath};
use itertools::Itertools;
use log::{debug, info, trace, warn};
use marlu::{
    io::{ms::MeasurementSetWriter, VisWritable},
    precession::{precess_time, PrecessionInfo},
    time::gps_millis_to_epoch,
    RADec,
};
use prettytable::{cell, format as prettyformat, row, table};
use std::{env, ffi::OsString, fmt::Debug, ops::Range};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use birli::{
            flags::flag_jones_array_existing,
        };
        use aoflagger_sys::{cxx_aoflagger_new};
    }
}
use birli::{
    context_to_jones_array, correct_cable_lengths, correct_geometry, get_antenna_flags,
    init_flag_array,
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib::{CorrelatorContext, GeometricDelaysApplied},
        LatLngHeight,
    },
    write_flags,
};

// TODO: fix too_many_arguments
#[allow(clippy::too_many_arguments)]
pub fn show_param_info(
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
) {
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
        let epoch = gps_millis_to_epoch(gps_time_ms);
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
    let avg_timesteps = (timestep_range.len() as f64 / avg_time as f64).ceil() as usize;
    let avg_int_time_s = int_time_s * avg_time as f64;
    info!(
        "Output duration:      {:.3}s = {:3} * {:.3}s{}",
        avg_timesteps as f64 * avg_int_time_s,
        avg_timesteps,
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
    let out_channel_count = (coarse_chan_range.len() as f64 * fine_chans_per_coarse as f64
        / avg_freq as f64)
        .ceil() as usize;
    let avg_fine_chan_width_khz = fine_chan_width_khz * avg_freq as f64;
    info!(
        "Output Bandwidth:     {:.3}MHz = {:9} * {:.3}kHz{}",
        out_bandwidth_mhz,
        out_channel_count,
        avg_fine_chan_width_khz,
        if avg_freq != 1 {
            format!(" ({}x)", avg_freq)
        } else {
            "".into()
        }
    );

    let first_epoch = gps_millis_to_epoch(context.timesteps[0].gps_time_ms);
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
}

fn main_with_args<I, T>(args: I)
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
            arg!(--"time-sel" "[WIP] Timestep index range (inclusive) to select")
                .value_names(&["MIN", "MAX"])
                .required(false),
            arg!(--"time-chunk" <STEPS> "[WIP] Process observation in chunks of <STEPS> timesteps.")
                .required(false),
            arg!(--"max-memory" <GIBIBYTES> "[WIP] The maximum amount of memory to use")
                .help_heading("INPUT")
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

            // corrections
            arg!(--"phase-centre" "Override Phase centre from metafits (degrees)")
                .value_names(&["RA", "DEC"])
                .required(false),
            arg!(--"pointing-centre" "Use pointing instead phase centre")
                .conflicts_with("phase-centre"),
            arg!(--"no-cable-delay" "Do not perform cable length corrections"),
            arg!(--"no-geometric-delay" "Do not perform geometric corrections")
                .alias("no-geom"),
            arg!(--"emulate-cotter" "Use Cotter's array position, not MWAlib's"),
            arg!(--"dry-run" "Just print the summary and exit"),
            arg!(--"no-draw-progress" "do not show progress bars"),

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

    let metafits_path = matches
        .value_of("metafits")
        .expect("--metafits must be a valid path");
    let fits_paths: Vec<&str> = matches.values_of("fits_paths").expect("--").collect();

    let context =
        CorrelatorContext::new(&metafits_path, &fits_paths).expect("unable to get mwalib context");
    debug!("mwalib correlator context:\n{}", &context);
    let coarse_chan_range = get_coarse_chan_range(&context).unwrap();
    let mut coarse_chan_flags = get_coarse_chan_flags(&context);
    let sel_timestep_range = match matches.values_of("time-sel") {
        Some(mut values) => {
            if let (Some(from), Some(to)) = (values.next(), values.next()) {
                let from = from.parse::<usize>().expect("cannot parse --time-sel from");
                assert!(
                    from < context.num_timesteps,
                    "invalid --time-sel from {}. must be < num_timesteps ({})",
                    from,
                    context.num_timesteps
                );
                let to = to.parse::<usize>().expect("cannot parse --time-sel to");
                assert!(
                    to >= from && to < context.num_timesteps,
                    "invalid --time-sel from {} to {}, must be < num_timesteps ({})",
                    from,
                    to,
                    context.num_timesteps
                );
                from..to + 1
            } else {
                panic!("invalid --time-sel <from> <to>, two values must be provided");
            }
        }
        _ => get_timestep_range(&context).unwrap(),
    };

    let num_timesteps_per_chunk: usize = match (
        matches.value_of("time-chunk"),
        matches.value_of("max-memory"),
    ) {
        (Some(_), Some(_)) => {
            panic!("you can't use --time-chunk and --max-memory at the same time");
        }
        (Some(steps_str), None) => steps_str.parse().unwrap_or_else(|_| {
            panic!(
                "unable to parse --time-chunk \"{}\" as an unsigned integer",
                steps_str
            )
        }),
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
            let num_bytes_per_scan =
                context.num_timestep_coarse_chan_bytes * context.num_gpubox_files;
            let num_bytes_total = context.num_timesteps * num_bytes_per_scan;
            if max_memory_bytes < num_bytes_total as f64 {
                if max_memory_bytes < num_bytes_per_scan as f64 {
                    panic!(
                        "--max-memory ({}B) too small to fit a single timestep ({}GiB)",
                        max_memory_bytes, num_bytes_per_scan
                    );
                }
                (max_memory_bytes / num_bytes_per_scan as f64).floor() as usize
            } else {
                sel_timestep_range.len()
            }
        }
        _ => sel_timestep_range.len(),
    };

    let mut timestep_flags = get_timestep_flags(&context);
    let mut antenna_flags = get_antenna_flags(&context);
    let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
    let mut fine_chan_flags = vec![false; context.metafits_context.num_corr_fine_chans_per_coarse];

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

    let baseline_flags = get_baseline_flags(&context, &antenna_flags);

    // ///////// //
    // Show info //
    // ///////// //

    show_param_info(
        &context,
        array_pos,
        phase_centre,
        &coarse_chan_range,
        &sel_timestep_range,
        &baseline_idxs,
        &coarse_chan_flags,
        &fine_chan_flags,
        &timestep_flags,
        &antenna_flags,
        &baseline_flags,
        avg_time,
        avg_freq,
    );

    if matches.is_present("dry-run") {
        info!("Dry run. No files will be written.");
        return;
    }

    let draw_progress = !matches.is_present("no-draw-progress");

    for unimplemented_option in &[
        "flag-init",
        "flag-init-steps",
        "flag-end",
        "flag-end-steps",
        "flag-edge-width",
        "flag-edge-chans",
        "flag-fine-chans",
        "flag-dc",
        "no-flag-dc",
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
        "time-sel",
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

    let mut uvfits_writer = matches.value_of("uvfits-out").map(|uvfits_out| {
        UvfitsWriter::from_mwalib(
            uvfits_out,
            &context,
            &sel_timestep_range,
            &coarse_chan_range,
            &baseline_idxs,
            Some(array_pos),
            Some(phase_centre),
            avg_time,
            avg_freq,
        )
        .expect("couldn't initialise uvfits writer")
    });
    let mut ms_writer = matches.value_of("ms-out").map(|ms_out| {
        let writer = MeasurementSetWriter::new(ms_out, phase_centre, Some(array_pos));
        writer
            .initialize_from_mwalib(
                &context,
                &sel_timestep_range,
                &coarse_chan_range,
                &baseline_idxs,
                avg_time,
                avg_freq,
            )
            .unwrap();
        writer
    });

    for mut timestep_chunk in &sel_timestep_range.chunks(num_timesteps_per_chunk) {
        let chunk_first_timestep = timestep_chunk.next().unwrap();
        let chunk_last_timestep = timestep_chunk.last().unwrap_or(chunk_first_timestep);
        let chunk_timestep_range: Range<usize> = chunk_first_timestep..chunk_last_timestep + 1;
        let flag_array = init_flag_array(
            &context,
            &chunk_timestep_range,
            &coarse_chan_range,
            Some(&timestep_flags),
            Some(&coarse_chan_flags),
            Some(&fine_chan_flags),
            Some(&baseline_flags),
        );

        #[allow(unused_mut)]
        let (mut jones_array, mut flag_array) = context_to_jones_array(
            &context,
            &chunk_timestep_range,
            &coarse_chan_range,
            Some(flag_array),
            draw_progress,
        )
        .unwrap();

        // perform cable delays if user has not disabled it, and they haven't aleady beeen applied.

        let no_cable_delays = matches.is_present("no-cable-delay");
        let cable_delays_applied = context.metafits_context.cable_delays_applied;
        if !cable_delays_applied && !no_cable_delays {
            info!(
                "Applying cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
            correct_cable_lengths(&context, &mut jones_array, &coarse_chan_range, false);
        } else {
            info!(
                "Skipping cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
        }

        cfg_if! {
            if #[cfg(feature = "aoflagger")] {
                if !matches.is_present("no-rfi") {
                    let aoflagger = unsafe { cxx_aoflagger_new() };
                    let default_strategy_filename = aoflagger.FindStrategyFileMWA();
                    let strategy_filename = matches.value_of("aoflagger-strategy").unwrap_or(&default_strategy_filename);
                    info!("flagging with strategy {}", strategy_filename);
                    flag_array = flag_jones_array_existing(
                        &aoflagger,
                        strategy_filename,
                        &jones_array,
                        Some(flag_array),
                        true,
                        draw_progress,
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
                correct_geometry(
                    &context,
                    &mut jones_array,
                    &chunk_timestep_range,
                    &coarse_chan_range,
                    Some(array_pos),
                    Some(phase_centre),
                    false,
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
            write_flags(&context, &flag_array, flag_template, &coarse_chan_range)
                .expect("unable to write flags");
        }

        // perform averaging
        let num_pols = context.metafits_context.num_visibility_pols;
        let flag_array = expand_flag_array(flag_array.view(), num_pols);
        let weight_factor = get_weight_factor(&context);
        let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        // output uvfits
        if let Some(uvfits_writer) = uvfits_writer.as_mut() {
            uvfits_writer
                .write_vis_mwalib(
                    jones_array.view(),
                    weight_array.view(),
                    flag_array.view(),
                    &context,
                    &chunk_timestep_range,
                    &coarse_chan_range,
                    &baseline_idxs,
                    avg_time,
                    avg_freq,
                    draw_progress,
                )
                .expect("unable to write uvfits");
        }

        // output ms
        if let Some(ms_writer) = ms_writer.as_mut() {
            ms_writer
                .write_vis_mwalib(
                    jones_array.view(),
                    weight_array.view(),
                    flag_array.view(),
                    &context,
                    &chunk_timestep_range,
                    &coarse_chan_range,
                    &baseline_idxs,
                    avg_time,
                    avg_freq,
                    draw_progress,
                )
                .expect("unable to write ms");
        }
    }

    // Finalise the uvfits writer.
    if let Some(uvfits_writer) = uvfits_writer {
        uvfits_writer
            .write_ants_from_mwalib(&context.metafits_context)
            .expect("couldn't write antenna table to uvfits");
    }
}

fn main() {
    env_logger::try_init().unwrap_or(());
    trace!("start main");
    main_with_args(env::args());
    trace!("end main");
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

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {

    use super::main_with_args;
    use approx::abs_diff_eq;
    use birli::{io::mwaf::FlagFileSet, Complex};
    use csv::StringRecord;
    use fitsio::errors::check_status as fits_check_status;
    use float_cmp::{approx_eq, F32Margin, F64Margin};
    use itertools::izip;
    use lazy_static::lazy_static;
    use lexical::parse;
    use marlu::{
        fitsio, fitsio_sys,
        mwalib::{
            CorrelatorContext, _get_required_fits_key, _open_fits, _open_hdu, fits_open,
            fits_open_hdu, get_required_fits_key,
        },
        rubbl_casatables::{Table, TableOpenMode},
    };
    use regex::Regex;
    use std::{
        collections::{BTreeMap, HashSet},
        path::PathBuf,
    };
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

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--no-draw-progress",
            "-f",
            mwaf_path_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let mut birli_flag_file_set = FlagFileSet::open(
            mwaf_path_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        let mut cotter_flag_file_set = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        assert_flagsets_eq!(
            context,
            birli_flag_file_set,
            cotter_flag_file_set,
            gpubox_ids
        );
    }

    #[test]
    fn aoflagger_outputs_uvfits() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1247842824.uvfits");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        assert!(uvfits_path.exists());

        assert!(uvfits_path.metadata().unwrap().len() > 0);
    }

    fn get_1254670392_avg_paths() -> (&'static str, [&'static str; 24]) {
        let metafits_path = "tests/data/1254670392_avg/1254670392.fixed.metafits";
        let gpufits_paths = [
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox02_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox03_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox04_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox05_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox06_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox07_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox08_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox09_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox10_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox11_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox12_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox13_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox14_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox15_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox16_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox17_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox18_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox19_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox20_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox21_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox22_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox23_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox24_00.fits",
        ];
        (metafits_path, gpufits_paths)
    }

    lazy_static! {
        static ref COMPLEX_REGEX: Regex = Regex::new(format!(
                r"^(?P<only_real>{0})$|^(?P<only_imag>{0})j$|^\((?P<complex_real>{0})\+?(?P<complex_imag>{0})j\)$",
                r"-?[\d\.]+(e-?\d+)?"
            ).as_str()
        ).unwrap();
    }

    fn compare_uvfits_with_csv(
        uvfits_path: PathBuf,
        expected_csv_path: PathBuf,
        vis_margin: F32Margin,
    ) {
        // Check both files are present
        assert!(uvfits_path.exists());
        assert!(expected_csv_path.exists());
        // Check both files are not empty
        assert!(uvfits_path.metadata().unwrap().len() > 0);
        assert!(expected_csv_path.metadata().unwrap().len() > 0);

        // Parse our expected CSV header
        // let expected_file = File::open(expected_csv_path).unwrap();
        let mut expected_reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .flexible(true)
            .trim(csv::Trim::All)
            .from_path(expected_csv_path)
            .unwrap();

        let headers = expected_reader.headers().unwrap();

        let keys = ["timestep", "baseline", "u", "v", "w", "pol", "type", "0"];

        let indices = parse_csv_headers(headers, &keys);

        let freq_start_header = indices.get("0").unwrap().to_owned();
        // let freq_end

        // Test the fits file has been correctly populated.
        let mut fptr = fits_open!(&uvfits_path).unwrap();
        let vis_hdu = fits_open_hdu!(&mut fptr, 0).unwrap();

        let pcount: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "PCOUNT").unwrap();
        let pzeros: Vec<f64> = (0..5)
            .map(|p_idx| {
                get_required_fits_key!(&mut fptr, &vis_hdu, format!("PZERO{}", p_idx + 1).as_str())
                    .unwrap()
            })
            .collect();
        let floats_per_pol: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS2").unwrap();
        let num_pols: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS3").unwrap();
        let num_fine_freq_chans: usize =
            get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS4").unwrap();
        let floats_per_complex = 2;

        let vis_len = num_fine_freq_chans * num_pols * floats_per_pol;
        assert!(vis_len > 0);

        let mut status = 0;
        let mut row_idx = 0;
        let mut obs_vis: Vec<f32> = vec![0.0; vis_len];
        let mut obs_group_params: Vec<f64> = vec![0.0; pcount];

        let pol_order = vec!["xx", "yy", "xy", "yx"];
        assert_eq!(num_pols, pol_order.len());

        let time_resolution = 1. / 1_000_000.;
        let mut times_seen = HashSet::<u64>::new();

        for record in expected_reader.records().filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(err) => panic!("{:?}", err),
        }) {
            let exp_group_params = ["u", "v", "w", "baseline", "timestep"]
                .iter()
                .map(|key| {
                    let value = &record[indices[&key.to_string()]];
                    value
                        .parse::<f64>()
                        .unwrap_or_else(|_| panic!("unable to parse {} -> {}", key, value))
                })
                .collect::<Vec<_>>();

            // Skip baseline(0,0)
            if exp_group_params[3] as i32 == 257 {
                continue;
            }

            let rec_type = record.get(indices[&String::from("type")]).unwrap();
            let pol = record.get(indices[&String::from("pol")]).unwrap();
            let pol_idx = pol_order.iter().position(|x| *x == pol).unwrap();

            let mut match_found = false;

            // iterate over rows in the uvfits file until we find an approximate match on timestep / baseline
            while row_idx < vis_len {
                unsafe {
                    // ffggpe = fits_read_grppar_flt
                    fitsio_sys::ffggpd(
                        fptr.as_raw(),                 /* I - FITS file pointer                       */
                        1 + row_idx as i64, /* I - group to read (1 = 1st group)           */
                        1,                  /* I - first vector element to read (1 = 1st)  */
                        pcount as i64,      /* I - number of values to read                */
                        obs_group_params.as_mut_ptr(), /* O - array of values that are returned       */
                        &mut status, /* IO - error status                           */
                    );
                }
                fits_check_status(status).unwrap();

                for (value, pzero) in izip!(obs_group_params.iter_mut(), pzeros.iter()) {
                    *value += pzero
                }

                times_seen.insert((obs_group_params[4] / time_resolution).round() as u64);

                let time_match = approx_eq!(
                    f64,
                    exp_group_params[4],
                    obs_group_params[4],
                    F64Margin::default().epsilon(1e-5)
                );

                let baseline_match =
                    exp_group_params[3].round() as i32 == obs_group_params[3].round() as i32;

                if time_match && baseline_match {
                    match_found = true;

                    // Assert that the group params are equal
                    for (param_idx, (obs_group_param, exp_group_param)) in
                        izip!(obs_group_params.iter(), exp_group_params.iter()).enumerate()
                    {
                        assert!(
                            approx_eq!(
                                f64,
                                *obs_group_param,
                                *exp_group_param,
                                F64Margin::default().epsilon(1e-7)
                            ),
                            "cells don't match in param {}, row {}. {:?} != {:?}",
                            param_idx,
                            row_idx,
                            obs_group_params,
                            exp_group_params
                        );
                    }

                    unsafe {
                        // ffgpve = fits_read_img_flt
                        fitsio_sys::ffgpve(
                            fptr.as_raw(),        /* I - FITS file pointer                       */
                            1 + row_idx as i64,   /* I - group to read (1 = 1st group)           */
                            1,                    /* I - first vector element to read (1 = 1st)  */
                            obs_vis.len() as i64, /* I - number of values to read                */
                            0.0,                  /* I - value for undefined pixels              */
                            obs_vis.as_mut_ptr(), /* O - array of values that are returned       */
                            &mut 0,               /* O - set to 1 if any values are null; else 0 */
                            &mut status,          /* IO - error status                           */
                        );
                    };
                    fits_check_status(status).unwrap();

                    match rec_type {
                        "vis" => {
                            let exp_pol_vis: Vec<_> = record
                                .iter()
                                .skip(freq_start_header)
                                .flat_map(|cell| {
                                    let complex = parse_complex(cell);
                                    vec![complex.re, complex.im].into_iter()
                                })
                                .collect();

                            assert_eq!(
                                num_fine_freq_chans * num_pols * floats_per_complex,
                                exp_pol_vis.len() * num_pols
                            );

                            let obs_pol_vis: Vec<_> = obs_vis
                                .chunks(floats_per_pol * num_pols)
                                .flat_map(|chunk| {
                                    chunk.chunks(floats_per_pol).skip(pol_idx).take(1).flat_map(
                                        |complex_flag| {
                                            let conjugate = vec![complex_flag[0], -complex_flag[1]];
                                            conjugate
                                        },
                                    )
                                })
                                .collect();

                            for (vis_idx, (&obs_val, &exp_val)) in
                                izip!(obs_pol_vis.iter(), exp_pol_vis.iter()).enumerate()
                            {
                                assert!(
                                    approx_eq!(f32, obs_val, exp_val, vis_margin),
                                    "visibility cells don't match (obs {} != exp {}) in row {} (bl {} ts {}), pol {} ({}), vis index {}. \nobserved: {:?} != \nexpected: {:?}",
                                    obs_val,
                                    exp_val,
                                    row_idx,
                                    exp_group_params[3],
                                    exp_group_params[4],
                                    pol,
                                    pol_idx,
                                    vis_idx,
                                    &obs_pol_vis,
                                    &exp_pol_vis
                                );
                            }
                        }
                        "weight" => {
                            let exp_pol_weight: Vec<f32> = record
                                .iter()
                                .skip(freq_start_header)
                                .map(|cell| cell.parse().unwrap())
                                .collect();

                            assert_eq!(num_fine_freq_chans, exp_pol_weight.len());

                            let obs_pol_weight: Vec<_> = obs_vis
                                .chunks(floats_per_pol * num_pols)
                                .flat_map(|chunk| {
                                    chunk
                                        .chunks(floats_per_pol)
                                        .skip(pol_idx)
                                        .take(1)
                                        .map(|complex_flag| complex_flag[2])
                                })
                                .collect();

                            for (weight_idx, (&obs_val, &exp_val)) in
                                izip!(obs_pol_weight.iter(), exp_pol_weight.iter()).enumerate()
                            {
                                assert!(
                                    approx_eq!(f32, obs_val, exp_val, F32Margin::default()),
                                    "cells don't match (obs {} != exp {}) in row {} (bl {} ts {}), pol {} ({}), weight index {}. \nobserved: {:?} != \nexpected: {:?}",
                                    obs_val,
                                    exp_val,
                                    row_idx,
                                    exp_group_params[3],
                                    exp_group_params[4],
                                    pol,
                                    pol_idx,
                                    weight_idx,
                                    &obs_pol_weight,
                                    &exp_pol_weight
                                );
                            }
                        }
                        _ => {
                            panic!("unexpected record type {}", rec_type);
                        }
                    }

                    break;
                }

                row_idx += 1;
            }
            if !match_found {
                panic!(
                    "unable to find matching row for time={}, baseline={:?}, times_seen={:?}",
                    exp_group_params[4],
                    exp_group_params[3],
                    times_seen
                        .iter()
                        .map(|&x| (x as f64) * time_resolution)
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    fn compare_ms_with_csv(ms_path: PathBuf, expected_csv_path: PathBuf, vis_margin: F32Margin) {
        // Check both files are present
        assert!(ms_path.exists());
        assert!(expected_csv_path.exists());
        // Check both files are not empty
        assert!(ms_path.metadata().unwrap().len() > 0);
        assert!(expected_csv_path.metadata().unwrap().len() > 0);

        // Parse our expected CSV header
        // let expected_file = File::open(expected_csv_path).unwrap();
        let mut expected_reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .flexible(true)
            .trim(csv::Trim::All)
            .from_path(expected_csv_path)
            .unwrap();

        let headers = expected_reader.headers().unwrap();
        let keys = ["time", "ant1", "ant2", "u", "v", "w", "pol", "type", "0"];
        let indices = parse_csv_headers(headers, &keys);

        let freq_start_header = indices.get("0").unwrap().to_owned();

        // Test the ms file has been correctly populated.
        let mut main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        let num_rows = main_table.n_rows();
        // dbg!(num_rows);
        let data_tabledesc = main_table.get_col_desc("DATA").unwrap();
        let data_shape = data_tabledesc.shape().unwrap();
        // let num_freqs = data_shape[0] as usize;
        let num_pols = data_shape[1] as usize;

        let mut row_idx = 0;
        let mut mjds_seen = HashSet::<u64>::new();

        let pol_order = vec!["xx", "xy", "yx", "yy"];
        assert_eq!(num_pols, pol_order.len());

        for record in expected_reader.records().filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(err) => panic!("{:?}", err),
        }) {
            let exp_baseline: (usize, usize) = (
                record[indices["ant1"]].parse().unwrap(),
                record[indices["ant2"]].parse().unwrap(),
            );

            let exp_uvw: Vec<f64> = vec![
                record[indices["u"]].parse().unwrap(),
                record[indices["v"]].parse().unwrap(),
                record[indices["w"]].parse().unwrap(),
            ];

            let exp_mjd: f64 = record[indices["time"]].parse().unwrap();

            // Skip autos
            if exp_baseline.0 == exp_baseline.1 {
                continue;
            }

            let mut match_found = false;

            // iterate over rows in the ms file until we find an approximate match on timestep / baseline
            while row_idx < num_rows {
                // main_table.read_row(&mut row, row_idx).unwrap();

                let obs_mjd = main_table
                    .get_cell::<f64>("TIME_CENTROID", row_idx)
                    .unwrap();
                mjds_seen.insert((obs_mjd * 10.).round() as u64);
                // dbg!(exp_mjd, obs_mjd);
                let time_match =
                    approx_eq!(f64, exp_mjd, obs_mjd, F64Margin::default().epsilon(1e-5));

                let obs_baseline = (
                    main_table.get_cell::<i32>("ANTENNA1", row_idx).unwrap() as usize,
                    main_table.get_cell::<i32>("ANTENNA2", row_idx).unwrap() as usize,
                );

                let baseline_match = exp_baseline == obs_baseline as (usize, usize);

                if time_match && baseline_match {
                    match_found = true;

                    let obs_uvw = main_table.get_cell_as_vec::<f64>("UVW", row_idx).unwrap();
                    for (uvw_idx, (obs_uvw, exp_uvw)) in
                        izip!(obs_uvw.iter(), exp_uvw.iter()).enumerate()
                    {
                        assert!(
                            approx_eq!(f64, *obs_uvw, *exp_uvw, F64Margin::default().epsilon(1e-5)),
                            "cells don't match in UVW[{}], row {}. {:?} != {:?}",
                            uvw_idx,
                            row_idx,
                            obs_uvw,
                            exp_uvw
                        );
                    }

                    let rec_type = record.get(indices[&String::from("type")]).unwrap();
                    let pol = record.get(indices[&String::from("pol")]).unwrap();
                    let pol_idx = pol_order.iter().position(|x| *x == pol).unwrap();

                    match rec_type {
                        "vis" => {
                            let exp_pol_vis: Vec<Complex<f32>> = record
                                .iter()
                                .skip(freq_start_header)
                                .map(parse_complex)
                                .collect();

                            let obs_vis = main_table
                                .get_cell_as_vec::<Complex<f32>>("DATA", row_idx)
                                .unwrap();
                            let obs_pol_vis = obs_vis
                                .into_iter()
                                .skip(pol_idx)
                                .step_by(num_pols)
                                .collect::<Vec<_>>();
                            assert_eq!(obs_pol_vis.len(), exp_pol_vis.len());

                            for (vis_idx, (&obs_val, &exp_val)) in
                                izip!(obs_pol_vis.iter(), exp_pol_vis.iter()).enumerate()
                            {
                                assert!(
                                    abs_diff_eq!(obs_val, exp_val, epsilon = vis_margin.epsilon),
                                    "visibility arrays don't match (obs {} != exp {}) in row {} (bl {:?} ts {}), pol {} ({}), vis index {}. \nobserved: {:?} != \nexpected: {:?}",
                                    obs_val,
                                    exp_val,
                                    row_idx,
                                    exp_baseline,
                                    exp_mjd,
                                    pol,
                                    pol_idx,
                                    vis_idx,
                                    &obs_pol_vis,
                                    &exp_pol_vis
                                );
                            }
                        }
                        "weight" => {
                            let exp_pol_weight: Vec<f32> = record
                                .iter()
                                .skip(freq_start_header)
                                .map(|x| x.parse::<f32>().unwrap())
                                .collect();

                            let obs_weight = main_table
                                .get_cell_as_vec::<f32>("WEIGHT_SPECTRUM", row_idx)
                                .unwrap();

                            let obs_pol_weight = obs_weight
                                .into_iter()
                                .skip(pol_idx)
                                .step_by(num_pols)
                                .collect::<Vec<_>>();

                            assert_eq!(obs_pol_weight.len(), exp_pol_weight.len());

                            for (weight_idx, (&obs_val, &exp_val)) in
                                izip!(obs_pol_weight.iter(), exp_pol_weight.iter()).enumerate()
                            {
                                assert!(
                                    abs_diff_eq!(obs_val, exp_val, epsilon = vis_margin.epsilon),
                                    "weight arrays don't match (obs {} != exp {}) in row {} (bl {:?} ts {}), pol {} ({}), weight index {}. \nobserved: {:?} != \nexpected: {:?}",
                                    obs_val,
                                    exp_val,
                                    row_idx,
                                    exp_baseline,
                                    exp_mjd,
                                    pol,
                                    pol_idx,
                                    weight_idx,
                                    &obs_pol_weight,
                                    &exp_pol_weight
                                );
                            }
                        }
                        "flag" => {
                            let exp_pol_flag: Vec<bool> = record
                                .iter()
                                .skip(freq_start_header)
                                .map(|x| x.to_lowercase().parse::<bool>().unwrap())
                                .collect();

                            let obs_flag =
                                main_table.get_cell_as_vec::<bool>("FLAG", row_idx).unwrap();

                            let obs_pol_flag = obs_flag
                                .into_iter()
                                .skip(pol_idx)
                                .step_by(num_pols)
                                .collect::<Vec<_>>();

                            assert_eq!(obs_pol_flag.len(), exp_pol_flag.len());

                            for (flag_idx, (&obs_val, &exp_val)) in
                                izip!(obs_pol_flag.iter(), exp_pol_flag.iter()).enumerate()
                            {
                                assert!(
                                    obs_val == exp_val,
                                    "flag arrays don't match (obs {} != exp {}) in row {} (bl {:?} ts {}), pol {} ({}), flag index {}. \nobserved: {:?} != \nexpected: {:?}",
                                    obs_val,
                                    exp_val,
                                    row_idx,
                                    exp_baseline,
                                    exp_mjd,
                                    pol,
                                    pol_idx,
                                    flag_idx,
                                    &obs_pol_flag,
                                    &exp_pol_flag
                                );
                            }
                        }
                        _ => panic!("unexpected record type: {}", rec_type),
                    }

                    break;
                }

                row_idx += 1;
            }
            if !match_found {
                panic!(
                    "unable to find matching row for time={}, baseline={:?}, mjds_seen={:?}",
                    exp_mjd,
                    exp_baseline,
                    mjds_seen
                        .iter()
                        .map(|&x| (x as f64) / 10.)
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    fn parse_complex(cell: &str) -> Complex<f32> {
        let captures = COMPLEX_REGEX.captures(cell).unwrap();
        let (real, imag) = match (
            captures.name("complex_real"),
            captures.name("complex_imag"),
            captures.name("only_real"),
            captures.name("only_imag"),
        ) {
            (Some(real), Some(imag), _, _) => (
                parse::<f32, _>(real.as_str()).unwrap(),
                parse::<f32, _>(imag.as_str()).unwrap(),
            ),
            (None, None, Some(real), None) => (parse::<f32, _>(real.as_str()).unwrap(), 0.0),
            (None, None, None, Some(imag)) => (0.0, parse::<f32, _>(imag.as_str()).unwrap()),
            _ => panic!("can't parse complex {}", cell),
        };
        Complex::new(real, imag)
    }

    fn parse_csv_headers(headers: &StringRecord, keys: &[&str]) -> BTreeMap<String, usize> {
        let mut remaining_keys: HashSet<_> = keys.iter().map(|x| String::from(*x)).collect();
        let mut indices = BTreeMap::<String, usize>::new();

        for (idx, cell) in headers.iter().enumerate() {
            let mut remove: Option<String> = None;
            for key in remaining_keys.iter() {
                if cell == key {
                    indices.insert(String::from(cell), idx);
                    remove = Some(key.clone());
                    break;
                }
            }
            if let Some(key) = remove {
                remaining_keys.remove(&key);
            }
        }

        if !remaining_keys.is_empty() {
            panic!("not all keys found: {:?}", remaining_keys);
        }

        indices
    }

    #[test]
    fn test_1254670392_avg_uvfits_no_corrections() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");
        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.none.uvfits.csv");

        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-draw-progress",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default());
    }

    #[test]
    #[ignore = "slow"]
    fn test_1254670392_avg_uvfits_cable_only() {
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
            "--no-draw-progress",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.cable.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
        );
    }

    #[test]
    #[ignore = "slow"]
    fn test_1254670392_avg_uvfits_geom_only() {
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
            "--no-draw-progress",
            "--no-cable-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.geom.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
        );
    }

    #[test]
    fn test_1254670392_avg_uvfits_both() {
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
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path =
        //     PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.corrected.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
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
    fn test_1254670392_avg_uvfits_both_phase0() {
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
            "--phase-centre",
            "0.0",
            "0.0",
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-4),
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
    fn test_1254670392_avg_ms_phase_pointing() {
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
            "--pointing-centre",
            "-M",
            ms_path.to_str().unwrap(),
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
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
    fn test_1254670392_avg_ms_corrected() {
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
            "--no-draw-progress",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
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
    fn test_1254670392_avg_ms_none_avg_4s_160khz() {
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
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-res",
            "4",
            "--avg-freq-res",
            "160",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
        );
    }

    /// Same as above but with factors instead of resolution
    #[test]
    fn test_1254670392_avg_ms_none_avg_4s_160khz_factors() {
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
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor",
            "2",
            "--avg-freq-factor",
            "4",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
        );
    }

    /// Same as above but forcing chunks by using a small --max-memory
    #[test]
    fn test_1254670392_avg_ms_none_avg_4s_160khz_max_mem() {
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
            "--no-draw-progress",
            "--emulate-cotter",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--avg-time-factor", "2",
            "--avg-freq-factor", "4",
            "--time-sel", "0", "2",
            "--time-chunk", "2",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        let main_table = Table::open(&ms_path, TableOpenMode::Read).unwrap();
        assert_eq!(main_table.n_rows(), 2 * 8256);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
        );
    }

    // todo: test when time_chunk is not a multiple of avg_time
}
