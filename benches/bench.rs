use birli::{
    context_to_jones_array,
    flags::{expand_flag_array, flag_to_weight_array, get_weight_factor},
    io::write_uvfits,
    write_ms,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use glob::glob;
use marlu::mwalib::CorrelatorContext;
use ndarray::s;
use std::path::Path;
use std::{cmp::min, env, ops::Range};
use tempfile::tempdir;

fn get_test_dir() -> String {
    env::var("BIRLI_TEST_DIR").unwrap_or_else(|_| String::from("/mnt/data"))
}

const NUM_TIMESTEPS: usize = 10;

fn get_context_1196175296() -> CorrelatorContext {
    let test_dir = get_test_dir();
    let test_path = Path::new(&test_dir);
    let vis_path = test_path.join("1196175296_vis");
    let metafits_path = vis_path
        .join("1196175296.metafits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_glob = vis_path
        .join("1196175296_*gpubox*_00.fits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_files: Vec<String> = glob(gpufits_glob.as_str())
        .unwrap()
        .filter_map(Result::ok)
        .map(|path| path.to_str().unwrap().to_owned())
        .collect();
    CorrelatorContext::new(&metafits_path, &gpufits_files).unwrap()
}

/// Get the timestep_idxs, coarse_chan_idxs, and baseline_idxs for a given context.
fn get_indices(context: &CorrelatorContext) -> (Range<usize>, Range<usize>, Vec<usize>) {
    let first_common_timestep_idx = *(context.common_timestep_indices.first().unwrap());
    let last_timestep_idx = *(context.provided_timestep_indices.last().unwrap());
    // limit max number of timesteps
    let last_timestep_idx = min(first_common_timestep_idx + 10, last_timestep_idx);

    let img_timestep_range = first_common_timestep_idx..last_timestep_idx + 1;
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let img_coarse_chan_range =
        *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);
    // let mwalib_coarse_chan_range = 0..2;

    let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();
    (img_timestep_range, img_coarse_chan_range, img_baseline_idxs)
}

fn bench_uvfits_output_1196175296_none(crt: &mut Criterion) {
    let context = get_context_1196175296();

    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

    let (_, sel_coarse_chan_range, sel_baseline_idxs) = get_indices(&context);
    let sel_timestep_range = 0..NUM_TIMESTEPS;
    let (jones_array, flag_array) = context_to_jones_array(
        &context,
        &sel_timestep_range,
        &sel_coarse_chan_range,
        None,
        false,
    )
    .unwrap();

    let weight_factor = get_weight_factor(&context);
    let flag_array = expand_flag_array(flag_array.view(), 4);
    let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

    crt.bench_function(
        format!("uvfits_output - 1196175296 {} timesteps", NUM_TIMESTEPS).as_str(),
        |bch| {
            bch.iter(|| {
                write_uvfits(
                    black_box(uvfits_path.as_path()),
                    black_box(&context),
                    black_box(jones_array.slice(s![sel_timestep_range.clone(), .., ..])),
                    black_box(weight_array.slice(s![sel_timestep_range.clone(), .., .., ..])),
                    black_box(flag_array.slice(s![sel_timestep_range.clone(), .., .., ..])),
                    black_box(&sel_timestep_range),
                    black_box(&sel_coarse_chan_range),
                    black_box(&sel_baseline_idxs),
                    None,
                    None,
                    1,
                    1,
                    false,
                )
                .unwrap();
            });
        },
    );
}

fn bench_ms_output_1196175296_none(crt: &mut Criterion) {
    let context = get_context_1196175296();

    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.ms");

    let (_, sel_coarse_chan_range, sel_baseline_idxs) = get_indices(&context);
    let sel_timestep_range = 0..NUM_TIMESTEPS;
    let (jones_array, flag_array) = context_to_jones_array(
        &context,
        &sel_timestep_range,
        &sel_coarse_chan_range,
        None,
        false,
    )
    .unwrap();

    let weight_factor = get_weight_factor(&context);
    let flag_array = expand_flag_array(flag_array.view(), 4);
    let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

    crt.bench_function(
        format!("ms_output - 1196175296 {} timesteps", NUM_TIMESTEPS).as_str(),
        |bch| {
            bch.iter(|| {
                write_ms(
                    black_box(ms_path.as_path()),
                    black_box(&context),
                    black_box(jones_array.slice(s![sel_timestep_range.clone(), .., ..])),
                    black_box(weight_array.slice(s![sel_timestep_range.clone(), .., .., ..])),
                    black_box(flag_array.slice(s![sel_timestep_range.clone(), .., .., ..])),
                    black_box(&sel_timestep_range),
                    black_box(&sel_coarse_chan_range),
                    black_box(&sel_baseline_idxs),
                    None,
                    None,
                    1,
                    1,
                    false,
                )
                .unwrap();
            });
        },
    );
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_ms_output_1196175296_none,
        bench_uvfits_output_1196175296_none,
);
criterion_main!(benches);
