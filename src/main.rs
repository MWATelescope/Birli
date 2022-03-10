use birli::{
    cli::BirliContext,
    context_to_jones_array,
    flags::{add_dimension, flag_to_weight_array, get_weight_factor},
    io::{aocal::AOCalSols, WriteableVis},
    marlu::{
        io::{ms::MeasurementSetWriter, VisWritable},
        ndarray::s,
    },
    parse_args, with_increment_duration, Axis, BirliError, FlagFileSet, UvfitsWriter, VisSelection,
};
use itertools::Itertools;
use log::{info, trace};
use std::{collections::HashMap, env, ffi::OsString, fmt::Debug, time::Duration};

#[allow(clippy::field_reassign_with_default)]
fn main_with_args<I, T>(args: I)
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
    I: Debug,
{
    let birli_ctx = parse_args(args);
    if matches!(birli_ctx, Err(BirliError::DryRun {})) {
        info!("Dry run. No files will be written.");
        return;
    }

    let BirliContext {
        corr_ctx,
        prep_ctx,
        vis_sel,
        flag_ctx,
        io_ctx,
        avg_time,
        avg_freq,
        num_timesteps_per_chunk,
    } = birli_ctx.unwrap();

    // used to time large operations
    let mut durations = HashMap::<&str, Duration>::new();

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
    }

    let mut duration_sum = Duration::ZERO;
    for (name, duration) in durations {
        info!("{} duration: {:?}", name, duration);
        duration_sum += duration;
    }
    info!("total duration: {:?}", duration_sum);
}

fn main() {
    env_logger::init_from_env(
        env_logger::Env::default().filter_or(env_logger::DEFAULT_FILTER_ENV, "info"),
    );
    trace!("start main");
    main_with_args(env::args());
    trace!("end main");
}

#[cfg(test)]
mod tests {
    use std::env;

    use super::main_with_args;

    #[test]
    fn main_with_version_doesnt_crash() {
        main_with_args(&["birli", "--version"]);
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
            "--pfb-gains",
            "none",
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
    #[ignore = "chunks not supported for .mwaf"]
    fn aoflagger_outputs_flags_chunked() {
        let tmp_dir = tempdir().unwrap();
        let mwaf_path_template = tmp_dir.path().join("Flagfile%%.mwaf");

        let metafits_path = "tests/data/1247842824_flags/1247842824.metafits";
        let gpufits_paths =
            vec!["tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits"];

        #[rustfmt::skip]
        let mut args = vec![
            "birli",
            "-m",
            metafits_path,
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--time-chunk", "1",
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
                r"-?(nan|inf|[\d\.]+(e-?\d+)?)"
            ).as_str()
        ).unwrap();
    }

    fn compare_uvfits_with_csv(
        uvfits_path: PathBuf,
        expected_csv_path: PathBuf,
        vis_margin: F32Margin,
        ignore_weights: bool,
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
                            if ignore_weights {
                                break;
                            }
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

    fn compare_ms_with_csv(
        ms_path: PathBuf,
        expected_csv_path: PathBuf,
        vis_margin: F32Margin,
        ignore_weights: bool,
    ) {
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
                                if obs_val.is_nan() && exp_val.is_nan() {
                                    continue;
                                }
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
                            if ignore_weights {
                                break;
                            }
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
                            if ignore_weights {
                                break;
                            }
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
        let captures = COMPLEX_REGEX
            .captures(cell)
            .unwrap_or_else(|| panic!("invalid complex number: {}", cell));
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-cable-delay",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), false);
    }

    #[test]
    fn test_1254670392_avg_uvfits_none_chunked() {
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

        main_with_args(&args);

        compare_uvfits_with_csv(uvfits_path, expected_csv_path, F32Margin::default(), true);
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--no-geometric-delay",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths[..1]);

        main_with_args(&args);
        // let uvfits_path = PathBuf::from("/mnt/data/1254670392_vis/1254670392.birli.cable.uvfits");
        compare_uvfits_with_csv(
            uvfits_path,
            expected_csv_path,
            F32Margin::default().epsilon(5e-5),
            false,
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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
            false,
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
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

        main_with_args(&args);
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

        main_with_args(&args);
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(2e-4),
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
            "--no-digital-gains",
            "--no-draw-progress",
            "--pfb-gains",
            "none",
            "--emulate-cotter",
        ];
        args.extend_from_slice(&gpufits_paths);

        main_with_args(&args);

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

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
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

        main_with_args(&args);

        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-7),
            false,
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

        main_with_args(&args);

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
    fn test_1254670392_avg_ms_none_avg_4s_160khz_tiny_chunk() {
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

        main_with_args(&args);
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
    fn test_1254670392_avg_ms_none_norfi_cal() {
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

        main_with_args(&args);

        // ignoring weights because Cotter doesn't flag NaNs
        compare_ms_with_csv(ms_path, expected_csv_path, F32Margin::default(), true);
    }

    #[test]
    /// Handle when calibration solution is provided with 24 channels, but a subset of channels are provided
    fn test_1254670392_avg_ms_none_norfi_cal_partial() {
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

        main_with_args(&args);

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
    fn test_1254670392_avg_ms_none_norfi_nopfb() {
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

        main_with_args(&args);
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
    fn test_1254670392_avg_ms_none_norfi_nodigital() {
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

        main_with_args(&args);
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-3),
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
    fn test_1254670392_avg_ms_none_norfi_nodigital_pfb_cotter_40() {
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

        main_with_args(&args);
        compare_ms_with_csv(
            ms_path,
            expected_csv_path,
            F32Margin::default().epsilon(1e-2),
            false,
        );
    }
}
