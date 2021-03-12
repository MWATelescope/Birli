#[cxx::bridge]
#[allow(dead_code)]
pub mod ffi {

    unsafe extern "C++" {
        include!("birli/include/cxx_aoflagger.h");

        fn aoflagger_GetVersion(major: &mut i16, minor: &mut i16, subMinor: &mut i16);

        pub type CxxImageSet;
        pub type CxxFlagMask;
        pub type CxxAOFlagger;
        pub type CxxStrategy;
        unsafe fn cxx_aoflagger_new() -> UniquePtr<CxxAOFlagger>;

        // CxxImageSet methods
        fn Width(self: &CxxImageSet) -> usize;
        fn Height(self: &CxxImageSet) -> usize;
        fn ImageCount(self: &CxxImageSet) -> usize;
        fn HorizontalStride(self: &CxxImageSet) -> usize;
        // TODO: fix this
        #[allow(clippy::mut_from_ref)]
        fn ImageBuffer(self: &CxxImageSet, imgIndex: usize) -> &mut [f32];

        // CxxFlagMask methods
        fn Width(self: &CxxFlagMask) -> usize;
        fn Height(self: &CxxFlagMask) -> usize;
        fn HorizontalStride(self: &CxxFlagMask) -> usize;
        // TODO: fix this
        #[allow(clippy::mut_from_ref)]
        fn Buffer(self: &CxxFlagMask) -> &mut [bool];

        // CxxAOFlagger methods
        fn GetVersion(self: &CxxAOFlagger, major: &mut i16, minor: &mut i16, subMinor: &mut i16);
        // TODO: what if widthCapacity < width?
        unsafe fn MakeImageSet(
            self: &CxxAOFlagger,
            width: usize,
            height: usize,
            count: usize,
            initialValue: f32,
            widthCapacity: usize,
        ) -> UniquePtr<CxxImageSet>;
        unsafe fn MakeFlagMask(
            self: &CxxAOFlagger,
            width: usize,
            height: usize,
            initialValue: bool,
        ) -> UniquePtr<CxxFlagMask>;
        fn FindStrategyFileGeneric(self: &CxxAOFlagger, scenario: &String) -> String;
        fn FindStrategyFileMWA(self: &CxxAOFlagger) -> String;
        #[allow(clippy::ptr_arg)]
        fn LoadStrategyFile(self: &CxxAOFlagger, filename: &String) -> UniquePtr<CxxStrategy>;

        // CxxStrategy methods
        fn Run(self: &CxxStrategy, input: &CxxImageSet) -> UniquePtr<CxxFlagMask>;
        fn RunExisting(
            self: &CxxStrategy,
            input: &CxxImageSet,
            existingFlags: &CxxFlagMask,
        ) -> UniquePtr<CxxFlagMask>;
    }
}

#[cfg(test)]
mod tests {
    use super::ffi::{aoflagger_GetVersion, cxx_aoflagger_new};
    use std::mem::size_of_val;
    use std::os::raw::c_short;

    #[test]
    fn test_mem_layout() {
        let usize_test: usize = 0;
        assert_eq!(size_of_val(&usize_test), 8);
    }

    #[test]
    fn test_valid_aoflagger_version() {
        let mut major: c_short = -1;
        let mut minor: c_short = -1;
        let mut sub_minor: c_short = -1;

        #[allow(unused_unsafe)]
        unsafe {
            aoflagger_GetVersion(&mut major, &mut minor, &mut sub_minor);
        }
        assert!(major >= 3);
        assert!(minor >= 0);
        assert!(sub_minor >= 0);
    }

    #[test]
    fn test_valid_cxx_aoflagger_version() {
        let mut major: c_short = -1;
        let mut minor: c_short = -1;
        let mut sub_minor: c_short = -1;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            aoflagger.GetVersion(&mut major, &mut minor, &mut sub_minor);
        }
        assert!(major >= 3);
        assert!(minor >= 0);
        assert!(sub_minor >= 0);
    }

    #[test]
    fn test_valid_img_set_init() {
        let width = 2 as usize;
        let height = 3 as usize;
        let count = 4 as usize;
        let initial_val = 5 as f32;
        let width_cap = 6 as usize;
        let exp_stride = (((width_cap - 1) / 8) + 1) * 8;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let img_set = aoflagger.MakeImageSet(width, height, count, initial_val, width_cap);
            assert_eq!(img_set.Width(), width);
            assert_eq!(img_set.Height(), height);
            assert_eq!(img_set.ImageCount(), count);
            assert_eq!(img_set.HorizontalStride(), exp_stride);
            let fist_buffer = img_set.ImageBuffer(0);
            assert_eq!(fist_buffer[0], initial_val);
            // TODO: test for leaks, whether destructing is taking place?
        }
    }

    #[test]
    fn test_valid_img_set_rw() {
        let width = 2 as usize;
        let height = 3 as usize;
        let count = 4 as usize;
        let initial_val = 5 as f32;
        let width_cap = 6 as usize;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let img_set = aoflagger.MakeImageSet(width, height, count, initial_val, width_cap);
            let first_buffer_write = img_set.ImageBuffer(0);
            first_buffer_write[0] = 7 as f32;
            let first_buffer_read = img_set.ImageBuffer(0);
            assert_eq!(first_buffer_read[0], 7 as f32);
        }
    }

    #[test]
    fn test_valid_flag_mask_init() {
        let width = 21 as usize;
        let height = 22 as usize;
        let initial_val = false;
        let exp_stride = 8 * (width / 8 + 1);
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let flag_mask = aoflagger.MakeFlagMask(width, height, initial_val);
            assert_eq!(flag_mask.Width(), width);
            assert_eq!(flag_mask.Height(), height);
            assert_eq!(flag_mask.HorizontalStride(), exp_stride);
            let buffer = flag_mask.Buffer();
            assert_eq!(buffer[0], initial_val);
            // TODO: test for leaks, whether destructing is taking place?
        }
    }

    #[test]
    fn test_valid_flag_mask_rw() {
        let width = 21 as usize;
        let height = 22 as usize;
        let initial_val = false;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let flag_mask = aoflagger.MakeFlagMask(width, height, initial_val);
            let buffer_write = flag_mask.Buffer();
            buffer_write[0] = !initial_val;
            let buffer_read = flag_mask.Buffer();
            assert_eq!(buffer_read[0], !initial_val);
        }
    }

    #[test]
    fn test_find_strategy_file_generic() {
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file_default = aoflagger.FindStrategyFileGeneric(&String::from(""));
            assert!(str::ends_with(
                &strategy_file_default,
                "generic-default.lua"
            ));
            let strategy_file_minimal = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
            assert!(str::ends_with(
                &strategy_file_minimal,
                "generic-minimal.lua"
            ));
        }
    }

    #[test]
    fn test_find_strategy_file_mwa() {
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file = aoflagger.FindStrategyFileMWA();
            assert!(str::ends_with(&strategy_file, "mwa-default.lua"));
        }
    }

    #[test]
    fn test_load_strategy_file() {
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file = aoflagger.FindStrategyFileMWA();
            let strategy = aoflagger.LoadStrategyFile(&strategy_file);
            // TODO: better test
            assert_eq!(size_of_val(&strategy), 8);
        }
    }

    #[test]
    fn test_strategy_run() {
        let width = 5 as usize;
        let height = 6 as usize;
        let count = 4 as usize;
        let initial_val = 5 as f32;
        let width_cap = width as usize;
        let exp_stride = (((width - 1) / 4) + 1) * 4;
        let mut exp_flag_chunks = vec![vec![false; exp_stride]; height];
        let noise_x = 3;
        let noise_y = 4;
        let noise_z = 2;
        exp_flag_chunks[noise_y][noise_x] = true;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file_name = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
            let strategy = aoflagger.LoadStrategyFile(&strategy_file_name);
            let img_set = aoflagger.MakeImageSet(width, height, count, initial_val, width_cap);
            let img_buffer = img_set.ImageBuffer(noise_z);
            img_buffer[noise_y * exp_stride + noise_x] = 999 as f32;
            let flag_mask = strategy.Run(&img_set);
            let flag_stride = flag_mask.HorizontalStride();
            assert_eq!(flag_stride, exp_stride);
            let flag_buffer = flag_mask.Buffer();
            assert_eq!(
                &flag_buffer.chunks(exp_stride).collect::<Vec<_>>(),
                &exp_flag_chunks
            );
        }
    }

    #[test]
    /* TODO: is this an issue? https://github.com/MWATelescope/Birli/blob/dev/doc/aoflagger_strategy_binding_issue.md */
    #[ignore]
    fn test_strategy_run_existing() {
        let width = 5 as usize;
        let height = 6 as usize;
        let count = 4 as usize;
        let initial_val = 5 as f32;
        let width_cap = width as usize;
        let exp_stride = (((width - 1) / 4) + 1) * 4;
        let mut exp_flag_chunks = vec![vec![false; exp_stride]; height];
        let existing_x = 1;
        let existing_y = 2;
        exp_flag_chunks[existing_y][existing_x] = true;
        let noise_x = 3;
        let noise_y = 4;
        let noise_z = 2;
        exp_flag_chunks[noise_y][noise_x] = true;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file_name = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
            let strategy = aoflagger.LoadStrategyFile(&strategy_file_name);
            let img_set = aoflagger.MakeImageSet(width, height, count, initial_val, width_cap);
            let img_buffer = img_set.ImageBuffer(noise_z);
            img_buffer[noise_y * exp_stride + noise_x] = 999 as f32;
            let existing_flag_mask = aoflagger.MakeFlagMask(width, height, false);
            let existing_flag_buf = existing_flag_mask.Buffer();
            existing_flag_buf[existing_y * exp_stride + existing_x] = true;
            let flag_stride = existing_flag_mask.HorizontalStride();
            let flag_mask = strategy.RunExisting(&img_set, &existing_flag_mask);
            let flag_buffer = flag_mask.Buffer();
            assert_eq!(size_of_val(&flag_buffer), 16);

            assert_eq!(
                &flag_buffer.chunks(flag_stride).collect::<Vec<_>>(),
                &exp_flag_chunks
            );
        }
    }
}
