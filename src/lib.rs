mod cxx_aoflagger;
use cxx::UniquePtr;
// use cxx_aoflagger::ffi::{CxxAOFlagger, CxxFlagMask, CxxImageSet};
use cxx_aoflagger::ffi::{CxxAOFlagger, CxxImageSet};
// use mwalib::{Baseline, CorrelatorContext};
use mwalib::CorrelatorContext;
use std::collections::BTreeMap;

pub fn context_to_baseline_imgset_map(
    aoflagger: &CxxAOFlagger,
    context: &mut CorrelatorContext,
) -> BTreeMap<usize, UniquePtr<CxxImageSet>> {
    let coarse_chan_arr = context.coarse_chans.clone();
    let timestep_arr = context.timesteps.clone();

    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline =
        context.metafits_context.num_corr_fine_chans_per_coarse * floats_per_finechan;
    let height = context.num_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;
    let width = context.num_timesteps;
    let img_stride = (((width - 1) / 8) + 1) * 8;

    // let imgset = aoflagger.MakeImageSet(width, height, 8, 0 as f32, width);
    let mut baseline_imgset_map: BTreeMap<usize, UniquePtr<CxxImageSet>> = BTreeMap::new();
    for (baseline_idx, _) in context.metafits_context.baselines.iter().enumerate() {
        unsafe {
            baseline_imgset_map.insert(
                baseline_idx,
                aoflagger.MakeImageSet(width, height, 8, 0 as f32, width),
            );
        }
    }

    for (coarse_chan_idx, _) in coarse_chan_arr.iter().enumerate() {
        for (timestep_idx, _) in timestep_arr.iter().enumerate() {
            let img_buf = context
                .read_by_baseline(timestep_idx, coarse_chan_idx)
                .unwrap();
            for (baseline_idx, baseline_chunk) in img_buf.chunks(floats_per_baseline).enumerate() {
                for (fine_chan_idx, fine_chan_chunk) in
                    baseline_chunk.chunks(floats_per_finechan).enumerate()
                {
                    let x = timestep_idx;
                    let y = context.metafits_context.num_corr_fine_chans_per_coarse
                        * coarse_chan_idx
                        + fine_chan_idx;

                    let imgset = baseline_imgset_map.get(&baseline_idx).unwrap();
                    for (float_idx, float_val) in fine_chan_chunk.iter().enumerate() {
                        imgset.ImageBuffer(float_idx)[y * img_stride + x] = *float_val
                    }
                }
            }
        }
    }

    return baseline_imgset_map;
}

// pub fn flag_imgsets(
//     aoflagger: &CxxAOFlagger,
//     baseline_imgset_map: BTreeMap<u64, UniquePtr<CxxImageSet>>,
// ) {
//     unimplemented!();
// }

// pub fn write_flags(
//     context: &mut CorrelatorContext,
//     baseline_flagmask_map: BTreeMap<u64, UniquePtr<CxxFlagMask>>,
//     filename_template: String,
//     channel_ids: Vec<usize>,
// ) {
//     unimplemented!();
// }

#[cfg(test)]
mod tests {
    use super::context_to_baseline_imgset_map;
    use crate::cxx_aoflagger::ffi::cxx_aoflagger_new;
    use mwalib::CorrelatorContext;

    fn get_mwax_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
        let gpufits_paths = vec![
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    fn get_mwa_ord_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    #[test]
    fn test_context_to_baseline_imgset_map_mwax() {
        let mut context = get_mwax_context();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let baseline_imgset_map = context_to_baseline_imgset_map(&aoflagger, &mut context);

            let imgset0 = baseline_imgset_map.get(&0).unwrap();

            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 0], 0x410000 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 1], 0x410100 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 2], 0x410200 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 3], 0x410300 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 7], 0.0);

            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 0], 0x410008 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 1], 0x410108 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 2], 0x410208 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 3], 0x410308 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 7], 0.0);

            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 0], 0x410400 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 1], 0x410500 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 2], 0x410600 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 3], 0x410700 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 7], 0.0);

            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 0], 0x410408 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 1], 0x410508 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 2], 0x410608 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 3], 0x410708 as f32);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 7], 0.0);

            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 0], 0x410001 as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 1], 0x410101 as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 2], 0x410201 as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 3], 0x410301 as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 7], 0.0);

            /* ... */
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 0], 0x410007 as f32);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 1], 0x410107 as f32);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 2], 0x410207 as f32);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 3], 0x410307 as f32);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 7], 0.0);


            let imgset2 = baseline_imgset_map.get(&2).unwrap();

            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 0], 0x410020 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 1], 0x410120 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 2], 0x410220 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 3], 0x410320 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 7], 0.0);

            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 0], 0x410028 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 1], 0x410128 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 2], 0x410228 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 3], 0x410328 as f32);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 4], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 5], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 6], 0.0);
            assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 7], 0.0);
        }
    }

    #[test]
    fn test_context_to_baseline_imgset_map_mwa_ord() {
        let mut context = get_mwa_ord_context();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let baseline_imgset_map = context_to_baseline_imgset_map(&aoflagger, &mut context);

            let imgset0 = baseline_imgset_map.get(&0).unwrap();

            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 0], 0x10c5be as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 1], 0x14c5be as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 2], 0x18c5be as f32);
            assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 3], 0x1cc5be as f32);

            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 0], 0x10c5bf as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 1], 0x14c5bf as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 2], 0x18c5bf as f32);
            assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 3], 0x1cc5bf as f32);

            let imgset2 = baseline_imgset_map.get(&2).unwrap();

            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 0], 0x10c57e as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 1], 0x14c57e as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 2], 0x18c57e as f32);
            assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 3], 0x1cc57e as f32);

            assert_eq!(imgset2.ImageBuffer(1)[0 * img_stride + 0], -0x10c57f as f32);
            assert_eq!(imgset2.ImageBuffer(1)[0 * img_stride + 1], -0x14c57f as f32);
            assert_eq!(imgset2.ImageBuffer(1)[0 * img_stride + 2], -0x18c57f as f32);
            assert_eq!(imgset2.ImageBuffer(1)[0 * img_stride + 3], -0x1cc57f as f32);
        }
    }

    // #[test]
    // fn test_write_flags_mwax() {
    //     let mut context = get_mwa_ord_context();
    //     let baseline_flagmask_map = BTreeMap<u64, UniquePtr<CxxFlagMask>>::new();
    //     unsafe {
    //         let aoflagger = cxx_aoflagger_new();
    //         const flagm
    //         baseline_flagmask_map.insert(
    //             0, aoflagger
    //         )
    //     }
    // }
}
