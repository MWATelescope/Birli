#[cxx::bridge]
#[allow(dead_code)]
mod ffi {
    unsafe extern "C++" {
        include!("birli/include/cxx_aoflagger.h");

        fn aoflagger_GetVersion(major: &mut i16, minor: &mut i16, subMinor: &mut i16);

        type CxxImageSet;
        type CxxFlagMask;
        type CxxAOFlagger;
        unsafe fn cxx_aoflagger_new() -> UniquePtr<CxxAOFlagger>;

        // CxxImageSet methods
        fn Width(self: &CxxImageSet) -> usize;
        fn Height(self: &CxxImageSet) -> usize;
        fn ImageCount(self: &CxxImageSet) -> usize;
        fn HorizontalStride(self: &CxxImageSet) -> usize;
        // TODO: fix this
        #[allow(clippy::mut_from_ref)]
        fn ImageBuffer(self: &CxxImageSet, imageIndex: usize) -> &mut [f32];

        // CxxFlagMask methods
        fn Width(self: &CxxFlagMask) -> usize;
        fn Height(self: &CxxFlagMask) -> usize;
        fn HorizontalStride(self: &CxxFlagMask) -> usize;
        // TODO: fix this
        #[allow(clippy::mut_from_ref)]
        fn Buffer(self: &CxxFlagMask) -> &mut [u8];

        // CxxAOFlagger methods
        fn GetVersion(self: &CxxAOFlagger, major: &mut i16, minor: &mut i16, subMinor: &mut i16);
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
        fn FindStrategyFile(self: &CxxAOFlagger) -> String;
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
    fn test_valid_image_set_init() {
        let width = 2 as usize;
        let height = 3 as usize;
        let count = 4 as usize;
        let initial_value = 5 as f32;
        let width_capacity = 6 as usize;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let image_set =
                aoflagger.MakeImageSet(width, height, count, initial_value, width_capacity);
            assert_eq!(image_set.Width(), width);
            assert_eq!(image_set.Height(), height);
            assert_eq!(image_set.ImageCount(), count);
            assert_eq!(image_set.HorizontalStride(), 8);
            let fist_buffer = image_set.ImageBuffer(0);
            assert_eq!(fist_buffer[0], initial_value);
            // TODO: test for leaks, whether destructing is taking place?
        }
    }

    #[test]
    fn test_valid_image_set_rw() {
        let width = 2 as usize;
        let height = 3 as usize;
        let count = 4 as usize;
        let initial_value = 5 as f32;
        let width_capacity = 6 as usize;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let image_set =
                aoflagger.MakeImageSet(width, height, count, initial_value, width_capacity);
            let first_buffer_write = image_set.ImageBuffer(0);
            first_buffer_write[0] = 7 as f32;
            let first_buffer_read = image_set.ImageBuffer(0);
            assert_eq!(first_buffer_read[0], 7 as f32);
        }
    }

    #[test]
    fn test_valid_flag_mask_init() {
        let width = 21 as usize;
        let height = 22 as usize;
        let initial_value = false;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let flag_mask = aoflagger.MakeFlagMask(width, height, initial_value);
            assert_eq!(flag_mask.Width(), width);
            assert_eq!(flag_mask.Height(), height);
            assert_eq!(flag_mask.HorizontalStride(), 24);
            let buffer = flag_mask.Buffer();
            assert_eq!(buffer[0], 0 as u8);
            // TODO: test for leaks, whether destructing is taking place?
        }
    }

    #[test]
    fn test_valid_flag_mask_rw() {
        let width = 21 as usize;
        let height = 22 as usize;
        let initial_value = false;
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let flag_mask = aoflagger.MakeFlagMask(width, height, initial_value);
            let buffer_write = flag_mask.Buffer();
            buffer_write[0] = 4 as u8;
            let buffer_read = flag_mask.Buffer();
            assert_eq!(buffer_read[0], 4 as u8);
        }
    }

    #[test]
    fn test_valid_strategy_file() {
        unsafe {
            let aoflagger = cxx_aoflagger_new();
            let strategy_file = aoflagger.FindStrategyFile();
            assert!(str::ends_with(&strategy_file, "mwa-default.lua"));
        }
    }
}
