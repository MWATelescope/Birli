use clap::{crate_authors, crate_description, crate_name, crate_version, App, Arg, SubCommand};
use log::{debug, info};
// use log::{log_enabled, Level};
use std::{env, ffi::OsString, fmt::Debug};

use birli::{
    context_to_baseline_imgsets, correct_cable_lengths, cxx_aoflagger_new, flag_imgsets_existing,
    get_antenna_flags, get_aoflagger_version_string, get_flaggable_timesteps,
    init_baseline_flagmasks, write_flags,
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
                .takes_value(true)
                .required(true) // TODO: specify a default that works with mwa-ord and mwax
                .help("Sets the template used to name flag files. Percents are substituted for the zero-prefixed GPUBox ID, which can be up to 3 characters log. Similar to -o in Cotter. Example: FlagFile%%%.mwaf")
        )
        .arg(
            Arg::with_name("no-cable-delay")
                .long("no-cable-delay")
                .takes_value(false)
                .required(false)
                .number_of_values(0)
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
        let flag_template = aoflagger_matches.value_of("flag-template").unwrap();
        let fits_files: Vec<&str> = aoflagger_matches.values_of("fits-files").unwrap().collect();
        let context = match CorrelatorContext::new(&metafits_path, &fits_files) {
            Ok(context) => context,
            Err(err) => panic!("unable to get mwalib context: {}", err),
        };
        debug!("mwalib correlator context:\n{}", &context);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_timestep_idxs = match get_flaggable_timesteps(&context) {
            Ok(timestep_idxs) => timestep_idxs,
            Err(err) => panic!("unable to determine flaggable timesteps: {}", err),
        };

        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(get_antenna_flags(&context)),
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
            correct_cable_lengths(&context, &mut baseline_imgsets);
        }

        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        flag_imgsets_existing(
            &aoflagger,
            &strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();
        write_flags(&context, baseline_flagmasks, flag_template, &gpubox_ids);
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
    use birli::{flag_io::FlagFileSet, get_aoflagger_version_string};
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

        let num_baselines = context.metafits_context.num_baselines;
        let num_flags_per_row = context.metafits_context.num_corr_fine_chans_per_coarse;
        let num_common_timesteps = context.num_common_timesteps;
        let num_rows = num_common_timesteps * num_baselines;
        let num_flags_per_timestep = num_baselines * num_flags_per_row;

        assert!(num_baselines > 0);
        assert!(num_rows > 0);
        assert!(num_flags_per_row > 0);

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
        let birli_chan_header_flags_raw = birli_flag_file_set.read_chan_header_flags_raw().unwrap();

        let mut cotter_flag_file_set = FlagFileSet::open(
            "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf",
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();
        let cotter_chan_header_flags_raw =
            cotter_flag_file_set.read_chan_header_flags_raw().unwrap();

        for gpubox_id in gpubox_ids {
            let (birli_header, birli_flags) = birli_chan_header_flags_raw.get(&gpubox_id).unwrap();
            let (cotter_header, cotter_flags) =
                cotter_chan_header_flags_raw.get(&gpubox_id).unwrap();
            assert_eq!(birli_header.obs_id, cotter_header.obs_id);
            assert_eq!(birli_header.num_channels, cotter_header.num_channels);
            assert_eq!(birli_header.num_ants, cotter_header.num_ants);
            // assert_eq!(birli_header.num_common_timesteps, cotter_header.num_common_timesteps);
            assert_eq!(birli_header.num_timesteps, num_common_timesteps);
            assert_eq!(birli_header.num_pols, cotter_header.num_pols);
            assert_eq!(birli_header.gpubox_id, cotter_header.gpubox_id);
            assert_eq!(birli_header.bytes_per_row, cotter_header.bytes_per_row);
            // assert_eq!(birli_header.num_rows, cotter_header.num_rows);
            assert_eq!(birli_header.num_rows, num_rows);

            // assert_eq!(birli_flags.len(), cotter_flags.len());
            assert_eq!(
                birli_flags.len(),
                num_common_timesteps * num_baselines * num_flags_per_row
            );

            izip!(
                birli_flags.chunks(num_flags_per_timestep),
                cotter_flags.chunks(num_flags_per_timestep)
            ).enumerate().for_each(|(common_timestep_idx, (birli_timestep_chunk, cotter_timestep_chunk))| {
                izip!(
                    context.metafits_context.baselines.iter(),
                    birli_timestep_chunk.chunks(num_flags_per_row),
                    cotter_timestep_chunk.chunks(num_flags_per_row)
                ).enumerate().for_each(|(baseline_idx, (baseline, birli_baseline_chunk, cotter_baseline_chunk))| {
                    if baseline.ant1_index == baseline.ant2_index {
                        return
                    }

                    assert_eq!(
                        birli_baseline_chunk, cotter_baseline_chunk,
                        "flag chunks for common timestep {}, baseline {} (ants {}, {}) do not match! \nbirli:\n{:?}\ncotter:\n{:?}", 
                        common_timestep_idx, baseline_idx, baseline.ant1_index, baseline.ant2_index, birli_baseline_chunk, cotter_baseline_chunk
                    )
                });
            });
        }
    }
}
