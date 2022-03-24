use crate::{
    marlu::{
        fitsio, fitsio_sys,
        mwalib::{
            _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
            get_required_fits_key, CorrelatorContext,
        },
        rubbl_casatables::{Table, TableOpenMode},
    },
    Complex,
};
use approx::abs_diff_eq;
use csv::StringRecord;
use fitsio::errors::check_status as fits_check_status;
use float_cmp::{approx_eq, F32Margin, F64Margin};
use itertools::izip;
use lazy_static::lazy_static;
use lexical::parse;
use regex::Regex;
use std::{
    collections::{BTreeMap, HashSet},
    path::{Path, PathBuf},
};

#[macro_export]
macro_rules! compare_jones {
    ($a:expr, $b:expr) => {
        assert_abs_diff_eq!(TestJones::<f32>::from($a), TestJones::<f32>::from($b));
    };
}

pub const fn get_1254670392_avg_paths() -> (&'static str, [&'static str; 24]) {
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

#[allow(dead_code)]
pub(crate) fn get_1254670392_avg_context() -> CorrelatorContext {
    let (metafits_path, gpufits_paths) = get_1254670392_avg_paths();
    CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
}

pub fn get_mwax_context() -> CorrelatorContext {
    let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
    let gpufits_paths = vec![
        "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
        "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
        "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
        "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
    ];
    CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
}

pub fn get_mwa_ord_context() -> CorrelatorContext {
    let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
    let gpufits_paths = vec![
        "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
    ];
    CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
}

/// Get a dummy MWA Ord `corr_ctx` with multiple holes in the data
///
/// The gpubox (batch, hdu) tuples look like this:
/// - ts is according to [`marlu::mwalib::correlatorContext`]
///
/// |                   | ts=0   | 1      | 2      | 3      | 4      |
/// | ----------------- | ------ | ------ | ------ | ------ | ------ |
/// | gpubox=00         | (0, 0) | (0, 1) | .      | (1, 0) | .      |
/// | 01                | .      | (0, 0) | (0, 1) | (1, 0) | (1, 1) |
pub fn get_mwa_ord_dodgy_context() -> CorrelatorContext {
    let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
    let gpufits_paths = vec![
        "tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits",
        "tests/data/1196175296_mwa_ord/limited_1/1196175296_20171201145540_gpubox01_01.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
    ];
    CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
}

/// Get a dummy MWA Ord `corr_ctx` with no overlapping timesteps
///
/// The gpubox (batch, hdu) tuples look like this:
///
/// | gpubox \ timestep | 0      | 1      | 2      | 3      |
/// | ----------------- | ------ | ------ | ------ | ------ |
/// | 00                | (0, 0) | (0, 1) | .      | .      |
/// | 01                | .      | .      | (0, 1) | (1, 0) |
pub fn get_mwa_ord_no_overlap_context() -> CorrelatorContext {
    let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
    let gpufits_paths = vec![
        "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
        "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
    ];
    CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
}

/// Get a dummy MWA Ord `corr_ctx` with no timesteps
pub fn get_mwa_ord_no_timesteps_context() -> CorrelatorContext {
    // let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
    // let gpufits_paths =
    //     vec!["tests/data/1196175296_mwa_ord/empty/1196175296_20171201145440_gpubox01_00.fits"];
    let mut corr_ctx = get_mwa_ord_no_overlap_context();
    corr_ctx.provided_timestep_indices = vec![];
    corr_ctx
}

lazy_static! {
    static ref COMPLEX_REGEX: Regex = Regex::new(format!(
            r"^(?P<only_real>{0})$|^(?P<only_imag>{0})j$|^\((?P<complex_real>{0})\+?(?P<complex_imag>{0})j\)$",
            r"-?(nan|inf|[\d\.]+(e-?\d+)?)"
        ).as_str()
    ).unwrap();
}

pub fn compare_uvfits_with_csv(
    uvfits_path: &Path,
    expected_csv_path: PathBuf,
    vis_margin: F32Margin,
    ignore_weights: bool,
    ignore_missing_chans: bool,
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
    let num_fine_freq_chans: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS4").unwrap();
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
                let value = &record[indices[&(*key).to_string()]];
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
                    1 + row_idx as i64,            /* I - group to read (1 = 1st group)           */
                    1,                             /* I - first vector element to read (1 = 1st)  */
                    pcount as i64,                 /* I - number of values to read                */
                    obs_group_params.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut status,                   /* IO - error status                           */
                );
            }
            fits_check_status(status).unwrap();

            for (value, pzero) in izip!(obs_group_params.iter_mut(), pzeros.iter()) {
                *value += pzero;
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

                        if !ignore_missing_chans {
                            assert_eq!(
                                num_fine_freq_chans * num_pols * floats_per_complex,
                                exp_pol_vis.len() * num_pols
                            );
                        }

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

                        if !ignore_missing_chans {
                            assert_eq!(num_fine_freq_chans, exp_pol_weight.len());
                        }

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
        assert!(
            match_found,
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

// TODO: make this less shitty
#[allow(clippy::cognitive_complexity)]
pub fn compare_ms_with_csv(
    ms_path: &Path,
    expected_csv_path: PathBuf,
    vis_margin: F32Margin,
    ignore_weights: bool,
    ignore_missing_chans: bool,
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
            let time_match = approx_eq!(f64, exp_mjd, obs_mjd, F64Margin::default().epsilon(1e-5));

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

                        if !ignore_missing_chans {
                            assert_eq!(obs_pol_vis.len(), exp_pol_vis.len());
                        }

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

                        if !ignore_missing_chans {
                            assert_eq!(obs_pol_weight.len(), exp_pol_weight.len());
                        }

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

                        let obs_flag = main_table.get_cell_as_vec::<bool>("FLAG", row_idx).unwrap();

                        let obs_pol_flag = obs_flag
                            .into_iter()
                            .skip(pol_idx)
                            .step_by(num_pols)
                            .collect::<Vec<_>>();

                        if !ignore_missing_chans {
                            assert_eq!(obs_pol_flag.len(), exp_pol_flag.len());
                        }

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
        assert!(
            match_found,
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
        for key in &remaining_keys {
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

    assert!(
        remaining_keys.is_empty(),
        "not all keys found: {:?}",
        remaining_keys
    );

    indices
}
