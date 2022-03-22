use birli::{
    flags::{flag_to_weight_array, get_weight_factor},
    io::write_uvfits,
    write_ms, VisSelection,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use glob::glob;
use marlu::mwalib::CorrelatorContext;
use ndarray::s;
use std::env;
use std::path::Path;
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

fn bench_uvfits_output_1196175296_none(crt: &mut Criterion) {
    let corr_ctx = get_context_1196175296();

    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    vis_sel.timestep_range = 0..NUM_TIMESTEPS;

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel
        .allocate_jones(corr_ctx.metafits_context.num_corr_fine_chans_per_coarse)
        .unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

    let weight_factor = get_weight_factor(&corr_ctx);
    let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

    crt.bench_function(
        format!("uvfits_output - 1196175296 {} timesteps", NUM_TIMESTEPS).as_str(),
        |bch| {
            bch.iter(|| {
                write_uvfits(
                    black_box(uvfits_path.as_path()),
                    black_box(&corr_ctx),
                    black_box(jones_array.slice(s![vis_sel.timestep_range.clone(), .., ..])),
                    black_box(weight_array.slice(s![vis_sel.timestep_range.clone(), .., ..])),
                    black_box(&vis_sel.timestep_range),
                    black_box(&vis_sel.coarse_chan_range),
                    black_box(&vis_sel.baseline_idxs),
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
    let corr_ctx = get_context_1196175296();

    let tmp_dir = tempdir().unwrap();
    let ms_path = tmp_dir.path().join("1254670392.none.ms");
    let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    vis_sel.timestep_range = 0..NUM_TIMESTEPS;
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel
        .allocate_jones(corr_ctx.metafits_context.num_corr_fine_chans_per_coarse)
        .unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

    let weight_factor = get_weight_factor(&corr_ctx);
    let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

    crt.bench_function(
        format!("ms_output - 1196175296 {} timesteps", NUM_TIMESTEPS).as_str(),
        |bch| {
            bch.iter(|| {
                write_ms(
                    black_box(ms_path.as_path()),
                    black_box(&corr_ctx),
                    black_box(jones_array.slice(s![vis_sel.timestep_range.clone(), .., ..])),
                    black_box(weight_array.slice(s![vis_sel.timestep_range.clone(), .., ..])),
                    black_box(&vis_sel.timestep_range),
                    black_box(&vis_sel.coarse_chan_range),
                    black_box(&vis_sel.baseline_idxs),
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
