mod cxx_aoflagger;
use cxx::UniquePtr;
use cxx_aoflagger::ffi::{CxxAOFlagger, CxxImageSet};
use mwalib::CorrelatorContext;

pub unsafe fn context_to_imgset(
    aoflagger: &CxxAOFlagger,
    context: &mut CorrelatorContext,
    baseline: usize,
) -> UniquePtr<CxxImageSet> {
    let coarse_chan_arr = context.coarse_chans.clone();
    let timestep_arr = context.timesteps.clone();

    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline =
        context.metafits_context.num_corr_fine_chans_per_coarse * floats_per_finechan;
    let height = context.num_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;
    let width = context.num_timesteps;
    let img_stride = (((width - 1) / 8) + 1) * 8;

    let imgset = aoflagger.MakeImageSet(width, height, 8, 0 as f32, width);

    for (coarse_chan_idx, _) in coarse_chan_arr.iter().enumerate() {
        for (timestep_idx, _) in timestep_arr.iter().enumerate() {
            let img_buf = context
                .read_by_baseline(timestep_idx, coarse_chan_idx)
                .unwrap();
            for (baseline_idx, baseline_chunk) in img_buf.chunks(floats_per_baseline).enumerate() {
                if baseline_idx != baseline {
                    continue;
                }
                for (fine_chan_idx, fine_chan_chunk) in
                    baseline_chunk.chunks(floats_per_finechan).enumerate()
                {
                    let x = timestep_idx;
                    let y = context.metafits_context.num_corr_fine_chans_per_coarse
                        * coarse_chan_idx
                        + fine_chan_idx;
                    for (float_idx, float_val) in fine_chan_chunk.iter().enumerate() {
                        imgset.ImageBuffer(float_idx)[y * img_stride + x] = *float_val
                    }
                }
            }
        }
    }

    return imgset;
}

#[cfg(test)]
mod tests {
    use super::context_to_imgset;
    use crate::cxx_aoflagger::ffi::cxx_aoflagger_new;
    use mwalib::CorrelatorContext;

    #[test]
    fn test_context_to_imgset_mwax() {
        let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
        let gpufits_paths = vec![
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
        ];
        let mut context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let imgset = context_to_imgset(&aoflagger, &mut context, 0);

            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 0], 0x410000 as f32);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 1], 0x410100 as f32);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 2], 0x410200 as f32);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 3], 0x410300 as f32);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[0 * img_stride + 7], 0.0);

            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 0], 0x410008 as f32);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 1], 0x410108 as f32);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 2], 0x410208 as f32);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 3], 0x410308 as f32);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[1 * img_stride + 7], 0.0);

            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 0], 0x410400 as f32);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 1], 0x410500 as f32);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 2], 0x410600 as f32);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 3], 0x410700 as f32);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[2 * img_stride + 7], 0.0);

            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 0], 0x410408 as f32);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 1], 0x410508 as f32);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 2], 0x410608 as f32);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 3], 0x410708 as f32);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(0)[3 * img_stride + 7], 0.0);

            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 0], 0x410001 as f32);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 1], 0x410101 as f32);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 2], 0x410201 as f32);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 3], 0x410301 as f32);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(1)[0 * img_stride + 7], 0.0);

            /* ... */

            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 0], 0x410007 as f32);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 1], 0x410107 as f32);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 2], 0x410207 as f32);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 3], 0x410307 as f32);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 4], 0.0);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 5], 0.0);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 6], 0.0);
            assert_eq!(imgset.ImageBuffer(7)[0 * img_stride + 7], 0.0);
        }
    }

    #[test]
    fn test_context_to_imgset_mwa_ord() {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        let mut context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let imgset = context_to_imgset(&aoflagger, &mut context, 0);

            assert_eq!(imgset.ImageBuffer(0)[0], 0x10c5be as f32);
            assert_eq!(imgset.ImageBuffer(0)[1], 0x14c5be as f32);
            assert_eq!(imgset.ImageBuffer(0)[2], 0x18c5be as f32);
            assert_eq!(imgset.ImageBuffer(0)[3], 0x1cc5be as f32);

            assert_eq!(imgset.ImageBuffer(1)[0], 0x10c5bf as f32);
            assert_eq!(imgset.ImageBuffer(1)[1], 0x14c5bf as f32);
            assert_eq!(imgset.ImageBuffer(1)[2], 0x18c5bf as f32);
            assert_eq!(imgset.ImageBuffer(1)[3], 0x1cc5bf as f32);
        }
    }
}
