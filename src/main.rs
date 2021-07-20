use clap::{crate_authors, crate_description, crate_name, crate_version, App, Arg, SubCommand};
use log::{debug, info, trace};
// use log::{log_enabled, Level};
use std::{env, ffi::OsString, fmt::Debug, path::Path};

use birli::{
    context_to_baseline_imgsets, correct_cable_lengths, cxx_aoflagger_new, flag_imgsets_existing,
    get_antenna_flags, get_aoflagger_version_string, get_flaggable_timesteps,
    init_baseline_flagmasks,
    io::write_uvfits,
    write_flags,
};
// use birli::util::{dump_flagmask, dump_imgset};
use mwalib::CorrelatorContext;

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
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(antenna_flags),
        );

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
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
            &strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        if let Some(flag_template) = flag_template {
            write_flags(
                &context,
                &baseline_flagmasks,
                flag_template,
                &img_coarse_chan_idxs,
            )
            .unwrap();
        }

        if let Some(uvfits_out) = uvfits_out {
            let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
            write_uvfits(
                Path::new(uvfits_out),
                &context,
                &baseline_idxs,
                &baseline_imgsets,
                &baseline_flagmasks,
                &img_timestep_idxs,
                &img_coarse_chan_idxs,
            ).unwrap();
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
    use itertools::izip;
    use mwalib::CorrelatorContext;
    use std::env;
    use tempfile::tempdir;

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
        let filename_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        let mut args = vec![
            "birli",
            "aoflagger",
            "-m",
            metafits_path,
            "-f",
            filename_template.to_str().unwrap(),
        ];
        args.extend_from_slice(&gpufits_paths);
        dbg!(&args);

        main_with_args(&args);

        let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        assert!(!gpubox_ids.is_empty());

        let mut birli_flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
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

    // TODO: test uvfits output with / without cable-delay
}
