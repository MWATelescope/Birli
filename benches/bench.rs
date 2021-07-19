use birli::{context_to_baseline_imgsets, correct_cable_lengths, cxx_aoflagger_new};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use glob::glob;
use mwalib::CorrelatorContext;
use std::env;
use std::path::Path;

fn get_context_mwax_half_1247842824() -> CorrelatorContext {
    let test_dir = env::var("BIRLI_TEST_DIR").unwrap();
    let test_path = Path::new(&test_dir);
    let vis_path = test_path.join("1247842824_vis");
    let metafits_path = vis_path
        .join("1247842824.metafits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_glob = vis_path
        .join("1247842824_*gpubox*_00.fits")
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

fn get_context_ord_half_1196175296() -> CorrelatorContext {
    let test_dir = env::var("BIRLI_TEST_DIR").unwrap();
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

fn bench_context_to_baseline_imgsets_mwax_half_1247842824(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_mwax_half_1247842824();
    crt.bench_function(
        "context_to_baseline_imgsets - mwax_half_1247842824",
        |bch| {
            bch.iter(|| {
                context_to_baseline_imgsets(
                    black_box(&aoflagger),
                    black_box(&context),
                    &context.common_coarse_chan_indices.clone(),
                    &context.common_timestep_indices.clone(),
                    None,
                )
            })
        },
    );
}

fn bench_context_to_baseline_imgsets_ord_half_1196175296(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_ord_half_1196175296();
    crt.bench_function("context_to_baseline_imgsets - ord_half_1196175296", |bch| {
        bch.iter(|| {
            context_to_baseline_imgsets(
                black_box(&aoflagger),
                black_box(&context),
                &context.common_coarse_chan_indices.clone(),
                &context.common_timestep_indices.clone(),
                None,
            )
        })
    });
}

fn bench_correct_cable_lengths_mwax_half_1247842824(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_mwax_half_1247842824();
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let mut baseline_imgsets = context_to_baseline_imgsets(
        &aoflagger,
        &context,
        &img_coarse_chan_idxs,
        &context.common_timestep_indices.clone(),
        None,
    );
    crt.bench_function("correct_cable_lengths - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            correct_cable_lengths(
                black_box(&context),
                black_box(&mut baseline_imgsets),
                black_box(img_coarse_chan_idxs),
            )
        })
    });
}

fn bench_correct_cable_lengths_ord_half_1196175296(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_ord_half_1196175296();
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let mut baseline_imgsets = context_to_baseline_imgsets(
        &aoflagger,
        &context,
        &img_coarse_chan_idxs,
        &context.common_timestep_indices.clone(),
        None,
    );
    crt.bench_function("correct_cable_lengths - ord_half_1196175296", |bch| {
        bch.iter(|| {
            correct_cable_lengths(
                black_box(&context),
                black_box(&mut baseline_imgsets),
                black_box(img_coarse_chan_idxs),
            )
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_context_to_baseline_imgsets_mwax_half_1247842824,
        bench_context_to_baseline_imgsets_ord_half_1196175296,
        bench_correct_cable_lengths_mwax_half_1247842824,
        bench_correct_cable_lengths_ord_half_1196175296,
);
criterion_main!(benches);
