use birli::{cli::BirliContext, BirliError};
use log::{info, trace};
use std::{env, ffi::OsString, fmt::Debug, time::Duration};

#[allow(clippy::field_reassign_with_default)]
fn main_with_args<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
    I: Debug,
{
    let birli_ctx = match BirliContext::from_args(args) {
        Ok(birli_ctx) => birli_ctx,
        Err(BirliError::DryRun {}) => {
            info!("Dry run. No files will be written.");
            return 0;
        }
        Err(BirliError::ClapError(inner)) => inner.exit(),
        Err(e) => {
            eprintln!("error parsing args: {:?}", e);
            return 1;
        }
    };
    match birli_ctx.run() {
        Ok(durations) => {
            info!(
                "total duration: {:?}",
                durations
                    .into_iter()
                    .fold(Duration::ZERO, |duration_sum, (name, duration)| {
                        info!("{} duration: {:?}", name, duration);
                        duration_sum + duration
                    })
            );
            0
        }
        Err(e) => {
            eprintln!("preprocessing error: {}", e);
            1
        }
    }
}

fn main() {
    env_logger::init_from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    );
    trace!("start main");
    let retcode = main_with_args(env::args());
    trace!("end main");
    std::process::exit(retcode);
}

#[cfg(test)]
mod tests {
    use std::env;

    use tempfile::tempdir;

    use super::main_with_args;

    #[test]
    fn main_with_version_doesnt_crash() {
        assert_eq!(main_with_args(&["birli", "--version"]), 0);
    }

    #[test]
    fn main_with_dry_run_doesnt_crash() {
        #[rustfmt::skip]
        assert_eq!(
            main_with_args(&[
                "birli",
                "-m", "tests/data/1254670392_avg/1254670392.fixed.metafits",
                "--dry-run",
                "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits"
            ]),
            0
        );
    }

    #[test]
    fn main_with_bad_arg_returns_1() {
        #[rustfmt::skip]
        assert_ne!(
            main_with_args(&[
                "birli",
                "-m", "tests/data/1254670392_avg/1254670392.fixed.metafits",
                "--avg-time-factor", "0",
                "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits"
            ]),
            0
        );
    }

    #[test]
    fn main_succesful_writes_uvfits() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1247842824.uvfits");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

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
        ];
        args.extend_from_slice(&gpufits_paths);

        assert_eq!(main_with_args(&args), 0);

        assert!(uvfits_path.exists());
        assert!(uvfits_path.metadata().unwrap().len() > 0);
    }

    #[test]
    fn main_gracefully_handle_munted_cal_file() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path = tmp_dir.path().join("1247842824.uvfits");

        let metafits_path = "tests/data/1254670392_avg/1254670392.fixed.metafits";
        let gpufits_paths =
            ["tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits"];

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
            "--apply-di-cal", "tests/data/1254670392_avg/1254690096.munted.bin"
        ];
        args.extend_from_slice(&gpufits_paths);

        assert_ne!(main_with_args(&args), 0);
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
}
