use birli::{
    cli::BirliContext,
    get_durations,
    BirliError::{ClapError, DryRun},
};
use clap::ErrorKind::{DisplayHelp, DisplayVersion};
use log::{info, trace};
use marlu::{Complex, Jones};
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
        Err(DryRun {}) => {
            info!("Dry run. No files will be written.");
            return 0;
        }
        Err(ClapError(inner)) => {
            // Swallow broken pipe errors
            trace!("clap error: {:?}", inner.kind());
            let _ = inner.print();
            match inner.kind() {
                DisplayHelp | DisplayVersion => return 0,
                _ => return 1,
            }
        }
        Err(e) => {
            eprintln!("error parsing args: {e}");
            return 1;
        }
    };

    // Calculations for timing info
    let fine_chans_per_coarse = birli_ctx
        .corr_ctx
        .metafits_context
        .num_corr_fine_chans_per_coarse;
    let mem_selected_bytes = birli_ctx.vis_sel.estimate_bytes_best(fine_chans_per_coarse);
    let num_sel_timesteps = birli_ctx.vis_sel.timestep_range.len();
    let num_sel_chans = birli_ctx.vis_sel.coarse_chan_range.len() * fine_chans_per_coarse;
    let num_sel_baselines = birli_ctx.vis_sel.baseline_idxs.len();
    let num_sel_pols = birli_ctx.corr_ctx.metafits_context.num_visibility_pols;
    let mem_per_timestep_gib =
        mem_selected_bytes as f64 / num_sel_timesteps as f64 / 1024.0_f64.powi(3);
    let num_avg_timesteps =
        (birli_ctx.vis_sel.timestep_range.len() as f64 / birli_ctx.avg_time as f64).ceil() as usize;
    let num_avg_chans = (birli_ctx.vis_sel.coarse_chan_range.len() as f64
        * fine_chans_per_coarse as f64
        / birli_ctx.avg_freq as f64)
        .ceil() as usize;
    let avg_mem_per_timestep_bytes = num_avg_chans
        * num_sel_baselines
        * num_sel_pols
        * (std::mem::size_of::<Complex<f32>>()
            + std::mem::size_of::<f32>()
            + std::mem::size_of::<bool>());
    let avg_mem_per_timestep_gib = (avg_mem_per_timestep_bytes) as f64 / 1024.0_f64.powi(3);

    let result = match birli_ctx.channel_range_sel.ranges.len() {
        1 => birli_ctx.run(),
        n if n > 1 => birli_ctx.run_ranges(),
        _ => unreachable!(),
    };
    match result {
        Ok(_) => {
            info!(
                "total duration: {:?}",
                get_durations().into_iter().fold(
                    Duration::ZERO,
                    |duration_sum, (name, duration)| {
                        info!("{} duration: {:?}", name, duration);
                        duration_sum + duration
                    }
                )
            );
            let read_time = get_durations()
                .get("read")
                .expect("no data was read")
                .as_secs_f32();
            let read_rate_mibs = (mem_selected_bytes / 1024_usize.pow(2)) as f32 / read_time;
            let write_time = get_durations()
                .get("write")
                .expect("no data was written")
                .as_secs_f32();
            let write_rate_mibs = (avg_mem_per_timestep_bytes * num_avg_timesteps
                / 1024_usize.pow(2)) as f32
                / write_time;

            info!(
                "Estimated data read     = {:5}ts * {:6}ch * {:6}bl * ({}<Jones<f32>> + {}<f32> + {}<bool>) = {:7.02} GiB @ {:8.03} MiB/s",
                num_sel_timesteps,
                num_sel_chans,
                num_sel_baselines,
                std::mem::size_of::<Jones<f32>>(),
                std::mem::size_of::<f32>(),
                std::mem::size_of::<bool>(),
                mem_per_timestep_gib * num_sel_timesteps as f64,
                read_rate_mibs,
            );

            info!(
                "Estimated data written  = {:5}ts * {:6}ch * {:6}bl * {:1}pol * ({}<c32> + {}<f32> + {}<bool>)  = {:7.02} GiB @ {:8.03} MiB/s",
                num_avg_timesteps,
                num_avg_chans,
                num_sel_baselines,
                num_sel_pols,
                std::mem::size_of::<Complex<f32>>(),
                std::mem::size_of::<f32>(),
                std::mem::size_of::<bool>(),
                avg_mem_per_timestep_gib * num_avg_timesteps as f64,
                write_rate_mibs
            );
            0
        }
        // TODO(Dev): different return codes for different errors
        Err(e) => {
            eprintln!("preprocessing error: {e}");
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
    use approx::assert_abs_diff_eq;
    use birli::mwalib::{
        _open_fits, _open_hdu, fits_open, fits_open_hdu, get_required_fits_key, CorrelatorContext,
        _get_required_fits_key,
    };
    use tempfile::tempdir;

    use super::main_with_args;

    #[test]
    fn main_with_version_succeeds() {
        assert_eq!(main_with_args(["birli", "--version"]), 0);
    }

    #[test]
    fn main_with_help_succeeds() {
        assert_eq!(main_with_args(["birli", "--help"]), 0);
    }

    #[test]
    fn main_with_dry_run_doesnt_crash() {
        #[rustfmt::skip]
        assert_eq!(
            main_with_args([
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
            main_with_args([
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
    fn main_succesful_writes_picket_uvfits() {
        let tmp_dir = tempdir().unwrap();
        let uvfits_path_0 = tmp_dir.path().join("1119683928.uvfits");
        let uvfits_path_1 = tmp_dir.path().join("1119683928_ch63.uvfits");
        let uvfits_path_2 = tmp_dir.path().join("1119683928_ch69-70.uvfits");

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox02_00.fits",
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox03_00.fits",
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox04_00.fits",
        ];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m", metafits_path,
            "-u", uvfits_path_0.to_str().unwrap(),
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--sel-chan-ranges", "1,2-3"
        ];
        args.extend_from_slice(&gpufits_paths);

        assert_eq!(main_with_args(&args), 0);

        assert!(uvfits_path_1.exists());
        assert!(uvfits_path_1.metadata().unwrap().len() > 0);

        assert!(uvfits_path_2.exists());
        assert!(uvfits_path_2.metadata().unwrap().len() > 0);

        // check frequencies are correct.
        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        let cc1 = &corr_ctx.coarse_chans[1];
        let cc2 = &corr_ctx.coarse_chans[2];
        let cc3 = &corr_ctx.coarse_chans[3];

        let mut uvfits_fptr_1 = fits_open!(uvfits_path_1.as_path()).unwrap();
        let uvfits_1_hdu_0 = fits_open_hdu!(&mut uvfits_fptr_1, 0).unwrap();
        let result_center_freq_1: f64 =
            get_required_fits_key!(&mut uvfits_fptr_1, &uvfits_1_hdu_0, "CRVAL4").unwrap();

        let mut uvfits_fptr_2 = fits_open!(uvfits_path_2.as_path()).unwrap();
        let uvfits_2_hdu_0 = fits_open_hdu!(&mut uvfits_fptr_2, 0).unwrap();
        let result_center_freq_2: f64 =
            get_required_fits_key!(&mut uvfits_fptr_2, &uvfits_2_hdu_0, "CRVAL4").unwrap();

        // see: https://wiki.mwatelescope.org/display/MP/MWA+Fine+Channel+Centre+Frequencies
        let offset_40khz = 15_000.;
        let expected_center_freq_1 = cc1.chan_centre_hz as f64 + offset_40khz;
        let expected_center_freq_2 =
            (cc2.chan_centre_hz + cc3.chan_centre_hz) as f64 / 2. + offset_40khz;

        assert_abs_diff_eq!(result_center_freq_1, expected_center_freq_1);
        assert_abs_diff_eq!(result_center_freq_2, expected_center_freq_2);
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
    fn main_successful_with_metrics_output() {
        let tmp_dir = tempdir().unwrap();
        let metrics_path = tmp_dir.path().join("test_metrics.fits");

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "--metrics-out", metrics_path.to_str().unwrap(),
            "-m", metafits_path,
            "--sel-ants", "1", "2",
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
        ];
        args.extend_from_slice(&gpufits_paths);

        // This should succeed with the fixed array indexing
        assert_eq!(main_with_args(&args), 0);

        // Verify metrics file was created
        assert!(metrics_path.exists());
        assert!(metrics_path.metadata().unwrap().len() > 0);
    }

    #[test]
    fn main_successful_with_metrics_output_sequential_ants() {
        let tmp_dir = tempdir().unwrap();
        let metrics_path = tmp_dir.path().join("test_metrics_sequential.fits");

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "--metrics-out", metrics_path.to_str().unwrap(),
            "-m", metafits_path,
            "--sel-ants", "0", "1",
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains", "none",
            "--no-cable-delay",
            "--no-geometric-delay",
        ];
        args.extend_from_slice(&gpufits_paths);

        // This should also succeed with sequential antenna selection
        assert_eq!(main_with_args(&args), 0);

        // Verify metrics file was created
        assert!(metrics_path.exists());
        assert!(metrics_path.metadata().unwrap().len() > 0);
    }
}
