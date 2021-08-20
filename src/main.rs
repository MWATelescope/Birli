use clap::{crate_authors, crate_description, crate_name, crate_version, App, Arg, SubCommand};
use log::{debug, info, trace};
// use log::{log_enabled, Level};
use std::{env, ffi::OsString, fmt::Debug, path::Path};

use birli::{
    context_to_baseline_imgsets, correct_cable_lengths,
    corrections::correct_geometry,
    cxx_aoflagger_new, flag_imgsets_existing, get_antenna_flags, get_aoflagger_version_string,
    get_flaggable_timesteps, init_baseline_flagmasks,
    io::write_uvfits,
    mwa_rust_core::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib, LatLngHeight,
    },
    write_flags,
};
// use birli::util::{dump_flagmask, dump_imgset};
use mwalib::{CorrelatorContext, GeometricDelaysApplied};

fn main_with_args<I, T>(args: I)
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
    I: Debug,
{
    debug!("args:\n{:?}", &args);

    let aoflagger_version = get_aoflagger_version_string();
    let aoflagger_subcommand = SubCommand::with_name("aoflagger")
        .about("flag visibilities with aoFlagger")
        .version(aoflagger_version.as_str())
        .arg(
            Arg::with_name("metafits")
                .short("m")
                .takes_value(true)
                .required(true)
                .help("Sets the metafits file."),
        )
        .arg(
            Arg::with_name("fits-files")
                .required(true)
                .multiple(true)
                // .last(true)
        )
        .arg(
            Arg::with_name("flag-template")
                .short("f")
                .takes_value(true)// TODO: specify a default that works with mwa-ord and mwax
                .help("Sets the template used to name flag files. Percents are substituted for the zero-prefixed GPUBox ID, which can be up to 3 characters log. Similar to -o in Cotter. Example: FlagFile%%%.mwaf")
        )
        .arg(
            Arg::with_name("uvfits-out")
                .short("u")
                .takes_value(true)
                .help("Filename for uvfits output. Similar to -o in Cotter. Example: 1196175296.uvfits")
        )
        .arg(
            Arg::with_name("no-cable-delay")
                .long("no-cable-delay")
                .takes_value(false)
                .required(false)
                .help("Do not perform cable length corrections.")
        )
        .arg(
            Arg::with_name("no-geometric-delay")
                .long("no-geometric-delay")
                .takes_value(false)
                .required(false)
                .help("Do not perform geometric length corrections.")
        )
        .arg(
            Arg::with_name("emulate-cotter")
                .long("emulate-cotter")
                .takes_value(false)
                .required(false)
                .help("Use Cotter's value for array position instead of MWAlib for direct comparison with Cotter.")
        )
        // TODO: implement specify flag strategy
        // .arg(
        //     Arg::with_name("flag-strategy")
        //         .help("Set the strategy filename, e.g. /usr/local/share/aoflagger/strategies/generic-minimal.lua")
        // )
        ;
    let matches = App::new(crate_name!())
        .version(crate_version!())
        .author(crate_authors!())
        .about(crate_description!())
        .subcommand(aoflagger_subcommand)
        .get_matches_from(args);

    debug!("arg matches:\n{:?}", &matches);

    if let Some(aoflagger_matches) = matches.subcommand_matches("aoflagger") {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let metafits_path = aoflagger_matches.value_of("metafits").unwrap();
        let flag_template = aoflagger_matches.value_of("flag-template");
        let uvfits_out = aoflagger_matches.value_of("uvfits-out");
        let fits_files: Vec<&str> = aoflagger_matches.values_of("fits-files").unwrap().collect();
        let context = match CorrelatorContext::new(&metafits_path, &fits_files) {
            Ok(context) => context,
            Err(err) => panic!("unable to get mwalib context: {}", err),
        };
        debug!("mwalib correlator context:\n{}", &context);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        trace!("img_coarse_chan_idxs: {:?}", img_coarse_chan_idxs);
        let img_timestep_idxs = match get_flaggable_timesteps(&context) {
            Ok(timestep_idxs) => timestep_idxs,
            Err(err) => panic!("unable to determine flaggable timesteps: {}", err),
        };
        trace!("img_timestep_idxs: {:?}", img_timestep_idxs);
        let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let antenna_flags = get_antenna_flags(&context);
        // trace!("antenna_flags: {:?}", antenna_flags);
        trace!(
            "antenna_flags: {:?}",
            antenna_flags
                .iter()
                .enumerate()
                .filter_map(|(idx, &flag)| {
                    if flag {
                        Some(idx)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        );

        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(antenna_flags),
        );

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(&mut baseline_flagmasks),
        );

        // perform cable delays if user has not disabled it, and they haven't aleady beeen applied.

        let no_cable_delays = aoflagger_matches.is_present("no-cable-delay");
        let cable_delays_applied = context.metafits_context.cable_delays_applied;
        if !cable_delays_applied && !no_cable_delays {
            debug!(
                "Applying cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
            correct_cable_lengths(&context, &mut baseline_imgsets, img_coarse_chan_idxs);
        } else {
            debug!(
                "Skipping cable delays. applied: {}, desired: {}",
                cable_delays_applied, !no_cable_delays
            );
        }

        let strategy_filename = &aoflagger.FindStrategyFileMWA();
        debug!("flagging with strategy {}", strategy_filename);
        flag_imgsets_existing(
            &aoflagger,
            strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        let array_pos = if aoflagger_matches.is_present("emulate-cotter") {
            Some(LatLngHeight {
                longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
                latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
                height_metres: COTTER_MWA_HEIGHT_METRES,
            })
        } else {
            None
        };

        // perform geometric delaysq if user has not disabled it, and they haven't aleady beeen applied.
        let no_geometric_delays = aoflagger_matches.is_present("no-geometric-delay");
        let geometric_delays_applied = context.metafits_context.geometric_delays_applied;

        match (geometric_delays_applied, no_geometric_delays) {
            (GeometricDelaysApplied::No, false) => {
                debug!(
                    "Applying geometric delays. applied: {:?}, desired: {}",
                    geometric_delays_applied, !no_geometric_delays
                );
                correct_geometry(
                    &context,
                    &baseline_idxs,
                    &mut baseline_imgsets,
                    img_coarse_chan_idxs,
                    &img_timestep_idxs,
                    array_pos.clone(),
                );
            }
            (..) => {
                debug!(
                    "Skipping geometric delays. applied: {:?}, desired: {}",
                    geometric_delays_applied, !no_geometric_delays
                );
            }
        };

        // output flags

        if let Some(flag_template) = flag_template {
            write_flags(
                &context,
                &baseline_flagmasks,
                flag_template,
                img_coarse_chan_idxs,
            )
            .unwrap();
        }

        // output uvfits

        if let Some(uvfits_out) = uvfits_out {
            write_uvfits(
                Path::new(uvfits_out),
                &context,
                &baseline_idxs,
                &baseline_imgsets,
                &baseline_flagmasks,
                &img_timestep_idxs,
                img_coarse_chan_idxs,
                array_pos,
            )
            .unwrap();
        }
    }
}

fn main() {
    env_logger::try_init().unwrap_or(());
    info!("start main");
    main_with_args(env::args());
    info!("end main");
}

#[cfg(test)]
mod tests {
    use super::main_with_args;
    use birli::{get_aoflagger_version_string, io::mwaf::FlagFileSet};
    use fitsio::errors::check_status as fits_check_status;
    use float_cmp::{approx_eq, F32Margin, F64Margin};
    use itertools::izip;
    use mwa_rust_core::{fitsio, fitsio_sys, mwalib};
    use mwalib::{
        CorrelatorContext, _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
        get_required_fits_key,
    };
    use regex::Regex;
    use std::{
        collections::{BTreeMap, HashSet},
        env,
        path::PathBuf,
    };
    use tempfile::tempdir;

    use lexical::parse;

    #[test]
    fn main_with_version_doesnt_crash() {
        main_with_args(&["--version"]);
    }

    #[test]
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

    #[test]
    fn forked_aoflagger_with_version_prints_version() {
        let pkg_name = env!("CARGO_PKG_NAME");
        let aoflagger_version = get_aoflagger_version_string();
        assert_cli::Assert::main_binary()
            .with_args(&["aoflagger", "--version"])
            .succeeds()
            .stdout()
            .contains(format!("{}-aoflagger {}\n", pkg_name, aoflagger_version).as_str())
            .unwrap();
    }

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
            "aoflagger",
            "-m",
            metafits_path,
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
            "aoflagger",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-cable-delay",
            "--no-geometric-delay",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        assert!(uvfits_path.exists());

        assert!(uvfits_path.metadata().unwrap().len() > 0);
    }

    fn get_1254670392_avg_paths() -> (&'static str, [&'static str; 24]) {
        let metafits_path = "tests/data/1254670392_avg/1254670392.metafits";
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

        let mut remaining_keys: HashSet<_> =
            ["timestep", "baseline", "u", "v", "w", "pol", "type", "0"]
                .iter()
                .map(|x| String::from(*x))
                .collect();
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
        let mut obs_idx = 0;
        let mut obs_vis: Vec<f32> = vec![0.0; vis_len];
        let mut obs_group_params: Vec<f64> = vec![0.0; pcount];

        let float_regex = r"-?[\d\.]+(e-?\d+)?";
        //let complex_regex = r"^(?P<real>-?[\d\.]+)|(?P<imag>[\+-]?[\d\.]+j)|\(?(?P<real>-?[\d\.]+)?(?P<imag>[\+-]?[\d\.]+j)?\)?$").unwrap();

        let complex_regex =
            Regex::new(format!(
                r"^(?P<only_real>{0})$|^(?P<only_imag>{0})j$|^\((?P<complex_real>{0})\+?(?P<complex_imag>{0})j\)$", 
                float_regex
            ).as_str()
        ).unwrap();

        let pol_order = vec!["xx", "yy", "xy", "yx"];
        assert_eq!(num_pols, pol_order.len());

        for record in expected_reader.records().filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(err) => panic!("{:?}", err),
        }) {
            let exp_group_params = ["u", "v", "w", "baseline", "timestep"]
                .iter()
                .map(|key| {
                    record
                        .get(indices[&key.to_string()])
                        .unwrap()
                        .parse::<f64>()
                        .unwrap()
                })
                .collect::<Vec<_>>();

            // Skip baseline(0,0)
            if exp_group_params[3] as i32 == 257 {
                continue;
            }

            let rec_type = record.get(indices[&String::from("type")]).unwrap();

            if rec_type != "vis" {
                continue;
            }

            // iterate over rows in the uvfits file until we find an approximate match on timestep / baseline
            while obs_idx < vis_len {
                unsafe {
                    // ffggpe = fits_read_grppar_flt
                    fitsio_sys::ffggpd(
                        fptr.as_raw(),                 /* I - FITS file pointer                       */
                        1 + obs_idx as i64, /* I - group to read (1 = 1st group)           */
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

                let time_match = approx_eq!(
                    f64,
                    exp_group_params[4],
                    obs_group_params[4],
                    F64Margin::default().epsilon(1e-1)
                );

                let baseline_match = approx_eq!(
                    f64,
                    exp_group_params[3],
                    obs_group_params[3],
                    F64Margin::default().epsilon(1e-1)
                );

                if time_match && baseline_match {
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
                            obs_idx,
                            obs_group_params,
                            exp_group_params
                        );
                    }

                    let exp_pol_vis: Vec<_> = record
                        .iter()
                        .skip(freq_start_header)
                        .flat_map(|cell| {
                            let captures = complex_regex.captures(cell).unwrap();
                            // let complex_real = captures.name("complex_real");
                            // let complex_imag = captures.name("complex_imag");
                            // let only_real = captures.name("only_real");
                            // let only_imag = captures.name("only_imag");
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
                                (None, None, Some(real), None) => {
                                    (parse::<f32, _>(real.as_str()).unwrap(), 0.0)
                                }
                                (None, None, None, Some(imag)) => {
                                    (0.0, parse::<f32, _>(imag.as_str()).unwrap())
                                }
                                _ => panic!("can't parse complex {}", cell),
                            };
                            vec![real, imag].into_iter()
                        })
                        .collect();

                    assert_eq!(
                        num_fine_freq_chans * num_pols * floats_per_complex,
                        exp_pol_vis.len() * num_pols
                    );

                    unsafe {
                        // ffgpve = fits_read_img_flt
                        fitsio_sys::ffgpve(
                            fptr.as_raw(),        /* I - FITS file pointer                       */
                            1 + obs_idx as i64,   /* I - group to read (1 = 1st group)           */
                            1,                    /* I - first vector element to read (1 = 1st)  */
                            obs_vis.len() as i64, /* I - number of values to read                */
                            0.0,                  /* I - value for undefined pixels              */
                            obs_vis.as_mut_ptr(), /* O - array of values that are returned       */
                            &mut 0,               /* O - set to 1 if any values are null; else 0 */
                            &mut status,          /* IO - error status                           */
                        );
                    };
                    fits_check_status(status).unwrap();

                    let pol = record.get(indices[&String::from("pol")]).unwrap();
                    let pol_idx = pol_order.iter().position(|x| *x == pol).unwrap();

                    let obs_pol_vis: Vec<_> = obs_vis
                        .chunks(floats_per_pol * num_pols)
                        .flat_map(|chunk| {
                            chunk.chunks(floats_per_pol).skip(pol_idx).take(1).flat_map(
                                |complex_flag| {
                                    // complex_flag[0..2].iter()
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
                            "cells don't match (obs {} != exp {}) in row {} (bl {} ts {}), pol {} ({}), vis index {}. \nobserved: {:?} != \nexpected: {:?}",
                            obs_val,
                            exp_val,
                            obs_idx,
                            exp_group_params[3],
                            exp_group_params[4],
                            pol,
                            pol_idx,
                            vis_idx,
                            &obs_pol_vis,
                            &exp_pol_vis
                        );
                    }
                    break;
                }

                obs_idx += 1;
            }
        }
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
            "aoflagger",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];

        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.none.uvfits");
        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default());
    }

    #[test]
    fn test_1254670392_avg_uvfits_cable_only() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.cable.uvfits.csv");

        let mut args = vec![
            "birli",
            "aoflagger",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
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
    fn test_1254670392_avg_uvfits_geom_only() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1254670392.uvfits");

        let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();

        let expected_csv_path =
            PathBuf::from("tests/data/1254670392_avg/1254670392.cotter.geom.uvfits.csv");

        let mut args = vec![
            "birli",
            "aoflagger",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
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
            "aoflagger",
            "-m",
            metafits_path,
            "-u",
            uvfits_path.to_str().unwrap(),
            "--emulate-cotter",
        ];

        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);
        // let uvfits_path =
        //     PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.corrected.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
        );
    }
}
