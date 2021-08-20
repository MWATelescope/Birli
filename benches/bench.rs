use birli::{
    context_to_baseline_imgsets, context_to_jones_array, correct_cable_lengths, correct_geometry,
    cxx_aoflagger_new, init_baseline_flagmasks, io::write_uvfits,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use glob::glob;
use mwa_rust_core::mwalib;
use mwalib::CorrelatorContext;
use std::env;
use std::path::Path;
use tempfile::tempdir;

fn get_test_dir() -> String {
    env::var("BIRLI_TEST_DIR").unwrap_or_else(|_| String::from("/mnt/data"))
}

fn get_context_mwax_half_1247842824() -> CorrelatorContext {
    let test_dir = get_test_dir();
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

fn get_context_1254670392_avg() -> CorrelatorContext {
    let test_dir = get_test_dir();
    let test_path = Path::new(&test_dir);
    let vis_path = test_path.join("1254670392_vis");
    let metafits_path = vis_path
        .join("1254670392.metafits")
        .to_str()
        .unwrap()
        .to_owned();
    let gpufits_glob = vis_path
        .join("1254670392_*gpubox*_00.fits")
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

fn bench_context_to_jones_array_mwax_half_1247842824(crt: &mut Criterion) {
    let context = get_context_mwax_half_1247842824();
    let timestep_idxs = context.common_timestep_indices.clone();
    let timestep_range = timestep_idxs[0]..timestep_idxs[timestep_idxs.len() - 1] + 1;
    let coarse_chan_idxs = context.common_coarse_chan_indices.clone();
    let coarse_chan_range = coarse_chan_idxs[0]..coarse_chan_idxs[coarse_chan_idxs.len() - 1] + 1;
    let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();
    crt.bench_function("context_to_jones_array - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            context_to_jones_array(
                black_box(&context),
                black_box(&timestep_range),
                black_box(&coarse_chan_range),
                black_box(img_baseline_idxs.as_slice()),
            )
        })
    });
}

fn bench_context_to_jones_array_ord_half_1196175296(crt: &mut Criterion) {
    let context = get_context_ord_half_1196175296();
    let timestep_idxs = context.common_timestep_indices.clone();
    let timestep_range = timestep_idxs[0]..timestep_idxs[timestep_idxs.len() - 1] + 1;
    let coarse_chan_idxs = context.common_coarse_chan_indices.clone();
    let coarse_chan_range = coarse_chan_idxs[0]..coarse_chan_idxs[coarse_chan_idxs.len() - 1] + 1;
    let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();
    crt.bench_function("context_to_jones_array - ord_half_1196175296", |bch| {
        bch.iter(|| {
            context_to_jones_array(
                black_box(&context),
                black_box(&timestep_range),
                black_box(&coarse_chan_range),
                black_box(img_baseline_idxs.as_slice()),
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
        img_coarse_chan_idxs,
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
        img_coarse_chan_idxs,
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

fn bench_correct_geometry_mwax_half_1247842824(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_mwax_half_1247842824();
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let img_timestep_idxs = &context.common_timestep_indices;
    let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
    let mut baseline_imgsets = context_to_baseline_imgsets(
        &aoflagger,
        &context,
        img_coarse_chan_idxs,
        &context.common_timestep_indices.clone(),
        None,
    );
    crt.bench_function("correct_geometry - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            correct_geometry(
                black_box(&context),
                black_box(&baseline_idxs),
                black_box(&mut baseline_imgsets),
                black_box(img_coarse_chan_idxs),
                black_box(img_timestep_idxs),
                None,
            );
        })
    });
}

fn bench_correct_geometry_ord_half_1196175296(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_ord_half_1196175296();
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let img_timestep_idxs = &context.common_timestep_indices;
    let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
    let mut baseline_imgsets = context_to_baseline_imgsets(
        &aoflagger,
        &context,
        img_coarse_chan_idxs,
        &context.common_timestep_indices.clone(),
        None,
    );
    crt.bench_function("correct_geometry - ord_half_1196175296", |bch| {
        bch.iter(|| {
            correct_geometry(
                black_box(&context),
                black_box(&baseline_idxs),
                black_box(&mut baseline_imgsets),
                black_box(img_coarse_chan_idxs),
                black_box(img_timestep_idxs),
                None,
            );
        })
    });
}

fn bench_uvfits_output_1254670392_avg_none(crt: &mut Criterion) {
    let aoflagger = unsafe { cxx_aoflagger_new() };
    let context = get_context_1254670392_avg();
    let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
    let img_timestep_idxs = &context.common_timestep_indices;
    let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
    let baseline_imgsets = context_to_baseline_imgsets(
        &aoflagger,
        &context,
        img_coarse_chan_idxs,
        &context.common_timestep_indices.clone(),
        None,
    );

    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1254670392.none.uvfits");

    let baseline_flagmasks = init_baseline_flagmasks(
        &aoflagger,
        &context,
        img_coarse_chan_idxs,
        img_timestep_idxs,
        None,
    );

    crt.bench_function("uvfits_output - 1254670392_avg", |bch| {
        bch.iter(|| {
            write_uvfits(
                black_box(uvfits_path.as_path()),
                black_box(&context),
                black_box(&baseline_idxs),
                black_box(&baseline_imgsets),
                black_box(&baseline_flagmasks),
                black_box(img_timestep_idxs),
                black_box(img_coarse_chan_idxs),
                None,
            )
            .unwrap();
        });
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_context_to_baseline_imgsets_mwax_half_1247842824,
        bench_context_to_baseline_imgsets_ord_half_1196175296,
        bench_context_to_jones_array_mwax_half_1247842824,
        bench_context_to_jones_array_ord_half_1196175296,
        bench_correct_cable_lengths_mwax_half_1247842824,
        bench_correct_cable_lengths_ord_half_1196175296,
        bench_correct_geometry_mwax_half_1247842824,
        bench_correct_geometry_ord_half_1196175296,
        bench_uvfits_output_1254670392_avg_none,
);
criterion_main!(benches);
