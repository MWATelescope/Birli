use birli::{
    correct_cable_lengths, correct_geometry,
    flags::{flag_to_weight_array, get_weight_factor},
    io::write_uvfits,
    VisSelection,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use glob::glob;
use marlu::mwalib::CorrelatorContext;
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

fn bench_read_mwalib_mwax_half_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    crt.bench_function("VisSelection::read_mwalib - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            vis_sel
                .read_mwalib(
                    black_box(&corr_ctx),
                    black_box(jones_array.view_mut()),
                    black_box(flag_array.view_mut()),
                    false,
                )
                .unwrap()
        })
    });
}

fn bench_read_mwalib_ord_half_1196175296(crt: &mut Criterion) {
    let corr_ctx = get_context_ord_half_1196175296();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    crt.bench_function("VisSelection::read_mwalib - ord_half_1196175296", |bch| {
        bch.iter(|| {
            vis_sel
                .read_mwalib(
                    black_box(&corr_ctx),
                    black_box(jones_array.view_mut()),
                    black_box(flag_array.view_mut()),
                    false,
                )
                .unwrap()
        })
    });
}

fn bench_correct_cable_lengths_mwax_half_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

    crt.bench_function("correct_cable_lengths - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            correct_cable_lengths(
                black_box(&corr_ctx),
                black_box(&mut jones_array),
                black_box(&vis_sel.coarse_chan_range),
                false,
            )
        })
    });
}

fn bench_correct_cable_lengths_ord_half_1196175296(crt: &mut Criterion) {
    let corr_ctx = get_context_ord_half_1196175296();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();
    crt.bench_function("correct_cable_lengths - ord_half_1196175296", |bch| {
        bch.iter(|| {
            correct_cable_lengths(
                black_box(&corr_ctx),
                black_box(&mut jones_array),
                black_box(&vis_sel.coarse_chan_range),
                false,
            )
        })
    });
}

fn bench_correct_geometry_mwax_half_1247842824(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();
    crt.bench_function("correct_geometry - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            correct_geometry(
                black_box(&corr_ctx),
                black_box(&mut jones_array),
                black_box(&vis_sel.timestep_range),
                black_box(&vis_sel.coarse_chan_range),
                None,
                None,
                false,
            );
        })
    });
}

fn bench_correct_geometry_ord_half_1196175296(crt: &mut Criterion) {
    let corr_ctx = get_context_ord_half_1196175296();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();
    crt.bench_function("correct_geometry - ord_half_1196175296", |bch| {
        bch.iter(|| {
            correct_geometry(
                black_box(&corr_ctx),
                black_box(&mut jones_array),
                black_box(&vis_sel.timestep_range),
                black_box(&vis_sel.coarse_chan_range),
                None,
                None,
                false,
            );
        })
    });
}

fn bench_uvfits_output_ord_half_1196175296_none(crt: &mut Criterion) {
    let corr_ctx = get_context_ord_half_1196175296();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1196175296.none.uvfits");

    let weight_factor = get_weight_factor(&corr_ctx);
    let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

    crt.bench_function("uvfits_output - ord_half_1196175296", |bch| {
        bch.iter(|| {
            write_uvfits(
                black_box(uvfits_path.as_path()),
                black_box(&corr_ctx),
                black_box(jones_array.view()),
                black_box(weight_array.view()),
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
    });
}

fn bench_uvfits_output_mwax_half_1247842824_none(crt: &mut Criterion) {
    let corr_ctx = get_context_mwax_half_1247842824();
    let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
    vis_sel
        .read_mwalib(
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

    let tmp_dir = tempdir().unwrap();
    let uvfits_path = tmp_dir.path().join("1247842824.none.uvfits");

    let weight_factor = get_weight_factor(&corr_ctx);
    let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

    crt.bench_function("uvfits_output - mwax_half_1247842824", |bch| {
        bch.iter(|| {
            write_uvfits(
                black_box(uvfits_path.as_path()),
                black_box(&corr_ctx),
                black_box(jones_array.view()),
                black_box(weight_array.view()),
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
    });
}

criterion_group!(
    name = expensive_benches;
    config = Criterion::default().sample_size(10);
    targets =
        bench_read_mwalib_mwax_half_1247842824,
        bench_read_mwalib_ord_half_1196175296,
        bench_correct_cable_lengths_mwax_half_1247842824,
        bench_correct_cable_lengths_ord_half_1196175296,
        bench_correct_geometry_mwax_half_1247842824,
        bench_correct_geometry_ord_half_1196175296,
        bench_uvfits_output_ord_half_1196175296_none,
        bench_uvfits_output_mwax_half_1247842824_none
);
criterion_main!(expensive_benches);
