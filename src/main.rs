use clap::{crate_authors, crate_description, crate_name, crate_version, App, Arg, SubCommand};
use std::{env, ffi::OsString};

use birli::{
    context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets, get_aoflagger_version_string,
    write_flags,
};
use mwalib::CorrelatorContext;

fn main_with_args<I, T>(args: I)
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
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
        )
        .arg(
            Arg::with_name("flag-template")
                .short("f")
                .takes_value(true)
                .required(true) // TODO: specify a default that works with mwa-ord and mwax
                .help("Sets the template used to name flag files. Percents are substituted for the zero-prefixed GPUBox ID, which can be up to 3 characters log. Similar to -o in Cotter. Example: FlagFile%%%.mwaf")
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
    dbg!(&matches);

    if let Some(aoflagger_matches) = matches.subcommand_matches("aoflagger") {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let metafits_path = aoflagger_matches.value_of("metafits").unwrap();
        let flag_template = aoflagger_matches.value_of("flag-template").unwrap();
        let fits_files: Vec<&str> = aoflagger_matches.values_of("fits-files").unwrap().collect();
        let mut context = CorrelatorContext::new(&metafits_path, &fits_files).unwrap();
        let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &mut context);
        let strategy = aoflagger.LoadStrategyFile(&aoflagger.FindStrategyFileMWA());
        let baseline_flagmasks = flag_imgsets(strategy, baseline_imgsets);
        let gpubox_ids: Vec<usize> = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();
        write_flags(&context, baseline_flagmasks, flag_template, &gpubox_ids);
    }
}

#[cfg(not(tarpaulin_include))]
fn main() {
    main_with_args(env::args())
}

#[cfg(test)]
mod tests {
    use super::{get_aoflagger_version_string, main_with_args};
    use assert_cli;
    use birli::FlagFileSet;
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

        let gpubox_ids: Vec<usize> = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();

        let mut flag_file_set =
            FlagFileSet::open(&context, filename_template.to_str().unwrap(), &gpubox_ids).unwrap();
        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();
        let (chan1_header, chan1_flags_raw) = chan_header_flags_raw.get(&gpubox_ids[0]).unwrap();
        assert_eq!(chan1_header.gpubox_id, gpubox_ids[0]);
        let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, context.num_timesteps);
        assert_eq!(num_baselines, context.metafits_context.num_baselines);
        assert_eq!(chan1_header.num_channels, num_fine_chans_per_coarse);
        assert_eq!(
            chan1_flags_raw.len(),
            chan1_header.num_timesteps * num_baselines * chan1_header.num_channels
        );

        let tests = [
            (0, 0, 0, i8::from(true)),
            (0, 0, 1, i8::from(true)),
            (0, 1, 0, i8::from(false)),
            (0, 1, 1, i8::from(false)),
            (0, 2, 0, i8::from(false)),
            (0, 2, 1, i8::from(false)),
        ];
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests.iter() {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx * num_fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                &chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }
}
