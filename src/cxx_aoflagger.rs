#[cxx::bridge]
// These are necessary for binding to aoflagger
#[allow(clippy::ptr_arg)]
#[allow(clippy::clippy::too_many_arguments)]
pub mod ffi {
    unsafe extern "C++" {
        include!("birli/include/cxx_aoflagger.h");

        /// CXX Wrapper for [`aoflagger::ImageSet`], a set of time-frequency 'images' which
        /// together contain data for one correlated baseline.
        ///
        /// Please see the aoflagger documentation for more details.
        ///
        /// [`aoflagger::ImageSet`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1ImageSet.html
        pub type CxxImageSet;
        /// CXX Wrapper for [`aoflagger::FlagMask`], a two-dimensional mask of bool flags.
        ///
        /// Please see the aoflagger documentation for more details.
        ///
        /// [`aoflagger::FlagMask`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1FlagMask.html
        pub type CxxFlagMask;
        /// CXX Wrapper for [`aoflagger::AOFlagger`], the main class for access to the flagger functionality.
        ///
        /// Please see the aoflagger documentation for more details.
        ///
        /// [`aoflagger::AOFlagger`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1AOFlagger.html
        pub type CxxAOFlagger;
        /// CXX Wrapper for [`aoflagger::Strategy`], a flagging strategy definition.
        ///
        /// Please see the aoflagger documentation for more details.
        ///
        /// [`aoflagger::Strategy`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1Strategy.html
        pub type CxxStrategy;

        /// Create a new [`CxxAOFlagger`] instance
        unsafe fn cxx_aoflagger_new() -> UniquePtr<CxxAOFlagger>;

        // CxxImageSet methods
        /// Get the width (number of timesteps) of the [`CxxImageSet`]
        fn Width(self: &CxxImageSet) -> usize;
        /// Get the height (number of coarse × fine frequency channels) of the [`CxxImageSet`]
        fn Height(self: &CxxImageSet) -> usize;
        /// Get the count (number of polarizations × complex components) of the [`CxxImageSet`]
        fn ImageCount(self: &CxxImageSet) -> usize;
        /// Get the total number of floats in one row of the [`CxxImageSet`]
        ///
        /// Row might have been padded to allow for SSE instructions and other optimizations.
        /// Therefore, one should add the horizontal stride to a data pointer to get the float in
        /// the next row (channel).
        fn HorizontalStride(self: &CxxImageSet) -> usize;
        /// (Immutably) access the raw float buffer at `imgIndex` in the [`CxxImageSet`]
        fn ImageBuffer(self: &CxxImageSet, imgIndex: usize) -> &[f32];
        /// (Mutably) access the raw float buffer at `imgIndex` in the [`CxxImageSet`]
        fn ImageBufferMut(self: Pin<&mut CxxImageSet>, imgIndex: usize) -> &mut [f32];
        /// (Mutably) access the raw float buffer at `imgIndex` in the [`CxxImageSet`] without pins
        /// TODO: document safety
        #[allow(clippy::mut_from_ref)]
        unsafe fn ImageBufferMutUnsafe(self: &CxxImageSet, imgIndex: usize) -> &mut [f32];
        // unsafe fn AllImageBufferMuts(self: Pin<&mut CxxImageSet>) -> &[&mut [f32]];

        // CxxFlagMask methods
        /// Get the width (number of timesteps) of the [`CxxFlagMask`]
        fn Width(self: &CxxFlagMask) -> usize;
        /// Get the height (number of coarse × fine frequency channels) of the [`CxxFlagMask`]
        fn Height(self: &CxxFlagMask) -> usize;
        /// Get the total number of bools in one row of the [`CxxFlagMask`]
        ///
        /// Row might have been padded to allow for SSE instructions and other optimizations.
        /// Therefore, one should add the horizontal stride to a data pointer to get the float in
        /// the next row (channel).
        fn HorizontalStride(self: &CxxFlagMask) -> usize;
        /// (Immutably) access the raw bool buffer of the [`CxxFlagMask`]
        fn Buffer(self: &CxxFlagMask) -> &[bool];
        /// (Mutably) access the raw bool buffer of the [`CxxFlagMask`]
        fn BufferMut(self: Pin<&mut CxxFlagMask>) -> &mut [bool];

        // CxxAOFlagger methods
        /// Get the AOFlagger library version number separated in major, minor and subminor fields.
        ///
        /// # Examples
        ///
        /// ```rust
        /// use birli::cxx_aoflagger_new;
        /// use std::os::raw::c_short;
        ///
        /// let mut major: c_short = -1;
        /// let mut minor: c_short = -1;
        /// let mut sub_minor: c_short = -1;
        /// let aoflagger = unsafe { cxx_aoflagger_new() };
        /// aoflagger.GetVersion(&mut major, &mut minor, &mut sub_minor);
        /// assert!(major >= 3);
        /// assert!(minor >= 0);
        /// assert!(sub_minor >= 0);
        /// ```
        fn GetVersion(self: &CxxAOFlagger, major: &mut i16, minor: &mut i16, subMinor: &mut i16);
        /// Create a new [`CxxImageSet`] with specified dimensions and initial value.
        ///
        /// # Undefined Behavior
        ///
        /// TODO: what if widthCapacity < width?
        unsafe fn MakeImageSet(
            self: &CxxAOFlagger,
            width: usize,
            height: usize,
            count: usize,
            initialValue: f32,
            widthCapacity: usize,
        ) -> UniquePtr<CxxImageSet>;
        /// Create a new [`CxxFlagMask`] with specified dimensions and initial value.
        unsafe fn MakeFlagMask(
            self: &CxxAOFlagger,
            width: usize,
            height: usize,
            initialValue: bool,
        ) -> UniquePtr<CxxFlagMask>;
        /// Find a Lua strategy filename for a Generic telescope.
        fn FindStrategyFileGeneric(self: &CxxAOFlagger, scenario: &String) -> String;
        /// Find a Lua strategy filename for the MWA telescope.
        fn FindStrategyFileMWA(self: &CxxAOFlagger) -> String;
        /// Load a strategy from disk.
        fn LoadStrategyFile(self: &CxxAOFlagger, filename: &String) -> UniquePtr<CxxStrategy>;

        // CxxStrategy methods
        /// Run the flagging strategy on the given [`CxxImageSet`].
        fn Run(self: &CxxStrategy, input: &CxxImageSet) -> UniquePtr<CxxFlagMask>;
        /// Run the flagging strategy on the given [`CxxImageSet`] with existing flags..
        ///
        /// # Undefined Behavior
        ///
        /// Not sure if this actually works. See: [`aoflagger::Strategy.Run` Binding issue]
        ///
        /// [`aoflagger::Strategy.Run` Binding issue]: https://github.com/MWATelescope/Birli/blob/main/doc/aoflagger_strategy_binding_issue.md
        fn RunExisting(
            self: &CxxStrategy,
            input: &CxxImageSet,
            existingFlags: &CxxFlagMask,
        ) -> UniquePtr<CxxFlagMask>;
    }
}

// TODO:
// The C++ implementation of CxxAOFlagger, CxxImageSet and CxxFlagMask are thread safe.
// https://cxx.rs/extern-c++.html
unsafe impl Send for ffi::CxxImageSet {}
unsafe impl Sync for ffi::CxxImageSet {}
unsafe impl Send for ffi::CxxFlagMask {}
unsafe impl Sync for ffi::CxxFlagMask {}
unsafe impl Send for ffi::CxxAOFlagger {}
unsafe impl Sync for ffi::CxxAOFlagger {}

#[cfg(test)]
// TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
#[allow(clippy::float_cmp)]
mod tests {

    use super::ffi::cxx_aoflagger_new;
    use std::mem::size_of_val;

    #[test]
    fn test_mem_layout() {
        let usize_test: usize = 0;
        assert_eq!(size_of_val(&usize_test), 8);
    }

    #[test]
    fn test_valid_img_set_init() {
        let width = 2_usize;
        let height = 3_usize;
        let count = 4_usize;
        let initial_val = 5_f32;
        let width_cap = 6_usize;
        let exp_stride = (((width_cap - 1) / 8) + 1) * 8;
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let img_set =
            unsafe { aoflagger.MakeImageSet(width, height, count, initial_val, width_cap) };
        assert_eq!(img_set.Width(), width);
        assert_eq!(img_set.Height(), height);
        assert_eq!(img_set.ImageCount(), count);
        assert_eq!(img_set.HorizontalStride(), exp_stride);
        let fist_buffer = img_set.ImageBuffer(0);
        assert_eq!(fist_buffer[0], initial_val);
        // TODO: test for leaks, whether destructing is taking place?
    }

    #[test]
    fn test_valid_img_set_rw() {
        let width = 2_usize;
        let height = 3_usize;
        let count = 4_usize;
        let initial_val = 5_f32;
        let width_cap = 6_usize;
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut img_set =
            unsafe { aoflagger.MakeImageSet(width, height, count, initial_val, width_cap) };
        let first_buffer_write = img_set.pin_mut().ImageBufferMut(0);
        first_buffer_write[0] = 7_f32;
        let first_buffer_read = img_set.ImageBuffer(0);
        assert_eq!(first_buffer_read[0], 7_f32);
    }

    #[test]
    fn test_valid_flag_mask_init() {
        let width = 21_usize;
        let height = 22_usize;
        let initial_val = false;
        let exp_stride = 8 * (width / 8 + 1);
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let flag_mask = unsafe { aoflagger.MakeFlagMask(width, height, initial_val) };
        assert_eq!(flag_mask.Width(), width);
        assert_eq!(flag_mask.Height(), height);
        assert_eq!(flag_mask.HorizontalStride(), exp_stride);
        let buffer = flag_mask.Buffer();
        assert_eq!(buffer[0], initial_val);
    }

    #[test]
    fn test_valid_flag_mask_rw() {
        let width = 21_usize;
        let height = 22_usize;
        let initial_val = false;
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut flag_mask = unsafe { aoflagger.MakeFlagMask(width, height, initial_val) };
        let buffer_write = flag_mask.pin_mut().BufferMut();
        buffer_write[0] = !initial_val;
        let buffer_read = flag_mask.Buffer();
        assert_eq!(buffer_read[0], !initial_val);
    }

    #[test]
    fn test_find_strategy_file_generic() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
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

    #[test]
    fn test_find_strategy_file_mwa() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let strategy_file = aoflagger.FindStrategyFileMWA();
        assert!(str::ends_with(&strategy_file, "mwa-default.lua"));
    }

    #[test]
    fn test_load_strategy_file() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let strategy_file = aoflagger.FindStrategyFileMWA();
        let strategy = aoflagger.LoadStrategyFile(&strategy_file);
        // TODO: better test
        assert_eq!(size_of_val(&strategy), 8);
    }

    #[test]
    fn test_strategy_run() {
        let width = 5_usize;
        let height = 6_usize;
        let count = 4_usize;
        let initial_val = 5_f32;
        let width_cap = width as usize;
        let exp_stride = (((width - 1) / 4) + 1) * 4;
        let mut exp_flag_chunks = vec![vec![false; exp_stride]; height];
        let noise_x = 3;
        let noise_y = 4;
        let noise_z = 2;
        exp_flag_chunks[noise_y][noise_x] = true;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut img_set =
            unsafe { aoflagger.MakeImageSet(width, height, count, initial_val, width_cap) };
        let strategy_file_name = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
        let strategy = aoflagger.LoadStrategyFile(&strategy_file_name);
        let img_buffer = img_set.pin_mut().ImageBufferMut(noise_z);
        img_buffer[noise_y * exp_stride + noise_x] = 999_f32;
        let flag_mask = strategy.Run(&img_set);
        let flag_stride = flag_mask.HorizontalStride();
        assert_eq!(flag_stride, exp_stride);
        let flag_buffer = flag_mask.Buffer();
        assert_eq!(
            &flag_buffer.chunks(exp_stride).collect::<Vec<_>>(),
            &exp_flag_chunks
        );
    }

    #[test]
    /* TODO: is this an issue? https://github.com/MWATelescope/Birli/blob/dev/doc/aoflagger_strategy_binding_issue.md */
    #[ignore]
    fn test_strategy_run_existing() {
        let width = 5_usize;
        let height = 6_usize;
        let count = 4_usize;
        let initial_val = 5_f32;
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

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let strategy_file_name = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
        let strategy = aoflagger.LoadStrategyFile(&strategy_file_name);
        let mut img_set =
            unsafe { aoflagger.MakeImageSet(width, height, count, initial_val, width_cap) };
        let img_buffer = img_set.pin_mut().ImageBufferMut(noise_z);
        img_buffer[noise_y * exp_stride + noise_x] = 999_f32;
        let mut existing_flag_mask = unsafe { aoflagger.MakeFlagMask(width, height, false) };
        let existing_flag_buf = existing_flag_mask.pin_mut().BufferMut();
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
