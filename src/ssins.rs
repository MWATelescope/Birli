//! SSINS (Sky-Subtract Incoherent Noise Spectra)
//!
//! This module implements the SSINS algorithm, which is a method for
//! subtracting incoherent noise from the visibilities.
//!
//! The algorithm is described in:
//! - <https://arxiv.org/abs/1906.01093>
//! - <https://ssins.readthedocs.io/en/latest/>

use crate::marlu::{
    fitsio::{
        images::{ImageDescription, ImageType},
        FitsFile,
    },
    mwalib::CorrelatorContext,
    ndarray::s,
    ndarray::{Array2, Array3, ArrayView3},
    Jones, VisSelection,
};

// SSINS (Sky-Subtract Incoherent Noise Spectra)
pub(crate) struct SSINS {
    pub diff_mean_amp_tfp: Array3<f32>, // (times-1, frequencies, polarizations)
    pub flag_array: Array2<bool>,       // (times-1, frequencies)
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
}

impl SSINS {
    pub(crate) fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
    ) -> Self {
        // incoherently averaged (along baseline) difference between the visibilities in time.
        // product is a 3D array of shape (num_timesteps - 1, num_freqs, num_pols=4)
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let mut diff_mean_amp_tfp = Array3::<f32>::zeros((num_timesteps - 1, num_freqs, 4));

        let flag_array = Array2::<bool>::default((num_timesteps - 1, num_freqs));

        for t in 0..num_timesteps - 1 {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    let jones_diff = jones_array_tfb[[t, f, b]] - jones_array_tfb[[t + 1, f, b]];

                    diff_mean_amp_tfp[[t, f, 0]] += jones_diff[0].norm();
                    diff_mean_amp_tfp[[t, f, 1]] += jones_diff[1].norm();
                    diff_mean_amp_tfp[[t, f, 2]] += jones_diff[2].norm();
                    diff_mean_amp_tfp[[t, f, 3]] += jones_diff[3].norm();
                }
            }
        }

        let total_samples = (num_timesteps - 1) * num_baselines;
        diff_mean_amp_tfp /= total_samples as f32;

        let timesteps_nodiff = &corr_ctx.timesteps[chunk_vis_sel.timestep_range.clone()]
            .iter()
            .map(|ts| ts.gps_time_ms as f64 / 1000.0)
            .collect::<Vec<_>>();
        // Compute the average time between adjacent timesteps
        let timesteps_diff: Vec<f64> = timesteps_nodiff[1..]
            .iter()
            .zip(&timesteps_nodiff[..num_timesteps - 1])
            .map(|(a, b)| (a + b) / 2.0)
            .collect();
        let integration_time_s = timesteps_diff[1] - timesteps_diff[0];

        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(
            &chunk_vis_sel.coarse_chan_range.clone().collect::<Vec<_>>(),
        );
        let freq_width_hz = all_freqs_hz[1] - all_freqs_hz[0];

        Self {
            diff_mean_amp_tfp,
            flag_array,
            start_time_gps_s: timesteps_diff[0],
            integration_time_s,
            start_freq_hz: all_freqs_hz[0],
            channel_width_hz: freq_width_hz,
        }
    }

    #[cfg(feature = "aoflagger")]
    pub(crate) fn flag(&mut self, strategy_filename: Option<String>) {
        use aoflagger_sys::cxx_aoflagger_new;

        use crate::flags::{amps_tfp_to_imageset, flag_baseline_view_to_flagmask};

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let imgset = amps_tfp_to_imageset(&aoflagger, self.diff_mean_amp_tfp.view());
        let flag_strategy = if let Some(strategy) = strategy_filename {
            aoflagger.LoadStrategyFile(&strategy.to_string())
        } else {
            let strategy_filename = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
            aoflagger.LoadStrategyFile(&strategy_filename)
        };
        // This lets us pass in our mutable flag array view to something not expecting a mutable.
        let flagmask = flag_baseline_view_to_flagmask(&aoflagger, self.flag_array.view());
        let new_flagmask = flag_strategy.RunExisting(&imgset, &flagmask);
        let flag_buf = new_flagmask.Buffer();
        let stride = new_flagmask.HorizontalStride();
        for (img_timestep_idx, mut flag_timestep_view) in
            self.flag_array.outer_iter_mut().enumerate()
        {
            for (img_chan_idx, mut flag_singular_view) in
                flag_timestep_view.outer_iter_mut().enumerate()
            {
                flag_singular_view.fill(flag_buf[img_chan_idx * stride + img_timestep_idx]);
            }
        }
    }

    pub(crate) fn save_to_fits(
        &self,
        fptr: &mut FitsFile,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Write one image per polarization
        let pol_names = ["XX", "YY", "XY", "YX"];
        let (num_times, num_freqs, num_pols) = self.diff_mean_amp_tfp.dim();

        for pol_idx in 0..num_pols {
            let dim = [num_times, num_freqs];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("SSINS_POL={}", pol_names[pol_idx]);
            let hdu = fptr.create_image(&extname, &image_description)?;

            // Write the data for this polarization
            let pol_data: Vec<f32> = self
                .diff_mean_amp_tfp
                .slice(s![.., .., pol_idx])
                .iter()
                .cloned()
                .collect();

            hdu.write_image(fptr, &pol_data)?;

            // Basic image info
            hdu.write_key(fptr, "BUNIT", "arbitrary")?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;

            // Time axis info
            hdu.write_key(fptr, "CTYPE1", "FREQ")?;
            hdu.write_key(fptr, "CRVAL1", self.start_freq_hz as f64)?;
            hdu.write_key(fptr, "CDELT1", self.channel_width_hz as f64)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "Hz")?;

            // Frequency axis info
            hdu.write_key(fptr, "CTYPE2", "TIME")?;
            hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s as f64)?;
            hdu.write_key(fptr, "CDELT2", self.integration_time_s as f64 * 1000.0)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT2", "s")?;

            // SSINS-specific metadata
            hdu.write_key(fptr, "POL", pol_names[pol_idx])?;
            hdu.write_key(fptr, "OBJECT", "SSINS")?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "SSINS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;

            // Add description
            hdu.write_key(
                fptr,
                "COMMENT",
                "SSINS (Sky-Subtract Incoherent Noise Spectra)",
            )?;
            hdu.write_key(
                fptr,
                "COMMENT",
                "Incoherently averaged difference between visibilities in time",
            )?;
            hdu.write_key(
                fptr,
                "COMMENT",
                "One image per polarization (XX, YY, XY, YX)",
            )?;
        }

        // Write flag array
        let flag_dim = self.flag_array.dim();
        let flag_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[flag_dim.0, flag_dim.1],
        };
        let flag_extname = "SSINS_FLAGS";
        let hdu = fptr.create_image(flag_extname, &flag_image_description)?;
        hdu.write_image(
            fptr,
            &self
                .flag_array
                .iter()
                .cloned()
                .map(|b| if b { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        )?;

        // Basic image info
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;

        // Time axis info
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz as f64)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz as f64)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // Frequency axis info
        hdu.write_key(fptr, "CTYPE2", "TIME")?;
        hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s as f64)?;
        hdu.write_key(fptr, "CDELT2", self.integration_time_s as f64 * 1000.0)?;
        hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT2", "s")?;

        // SSINS-specific metadata
        hdu.write_key(fptr, "OBJECT", "SSINS")?;
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "SSINS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        Ok(())
    }
}

pub(crate) struct EAVILS {
    pub nodiff_mean_amp_tfp: Array3<f32>, // (times, frequencies, polarizations)
    pub flag_array: Array2<bool>,         // (times, frequencies)
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
}

impl EAVILS {
    pub(crate) fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let mut nodiff_mean_amp_tfp = Array3::<f32>::zeros((num_timesteps, num_freqs, 4));

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    let jones_nodiff = jones_array_tfb[[t, f, b]];
                    nodiff_mean_amp_tfp[[t, f, 0]] += jones_nodiff[0].norm();
                    nodiff_mean_amp_tfp[[t, f, 1]] += jones_nodiff[1].norm();
                    nodiff_mean_amp_tfp[[t, f, 2]] += jones_nodiff[2].norm();
                    nodiff_mean_amp_tfp[[t, f, 3]] += jones_nodiff[3].norm();
                }
            }
        }

        let total_samples = num_timesteps * num_baselines;
        nodiff_mean_amp_tfp /= total_samples as f32;

        let flag_array = Array2::<bool>::default((num_timesteps, num_freqs));

        let timesteps_nodiff = &corr_ctx.timesteps[chunk_vis_sel.timestep_range.clone()]
            .iter()
            .map(|ts| ts.gps_time_ms as f64 / 1000.0)
            .collect::<Vec<_>>();
        let integration_time_s = timesteps_nodiff[1] - timesteps_nodiff[0];

        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(
            &chunk_vis_sel.coarse_chan_range.clone().collect::<Vec<_>>(),
        );
        let freq_width_hz = all_freqs_hz[1] - all_freqs_hz[0];

        Self {
            nodiff_mean_amp_tfp,
            flag_array,
            start_time_gps_s: timesteps_nodiff[0],
            integration_time_s,
            start_freq_hz: all_freqs_hz[0],
            channel_width_hz: freq_width_hz,
        }
    }

    pub(crate) fn save_to_fits(
        &self,
        fptr: &mut FitsFile,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Write one image per polarization
        let pol_names = ["XX", "YY", "XY", "YX"];
        let (num_times, num_freqs, num_pols) = self.nodiff_mean_amp_tfp.dim();

        for pol_idx in 0..num_pols {
            let dim = [num_times, num_freqs];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("EAVILS_POL={}", pol_names[pol_idx]);
            let hdu = fptr.create_image(&extname, &image_description)?;

            // Write the data for this polarization
            let pol_data: Vec<f32> = self
                .nodiff_mean_amp_tfp
                .slice(s![.., .., pol_idx])
                .iter()
                .cloned()
                .collect();

            hdu.write_image(fptr, &pol_data)?;

            // Basic image info
            hdu.write_key(fptr, "BUNIT", "arbitrary")?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;

            // Time axis info
            hdu.write_key(fptr, "CTYPE1", "FREQ")?;
            hdu.write_key(fptr, "CRVAL1", self.start_freq_hz as f64)?;
            hdu.write_key(fptr, "CDELT1", self.channel_width_hz as f64)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "Hz")?;

            // Frequency axis info
            hdu.write_key(fptr, "CTYPE2", "TIME")?;
            hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s as f64)?;
            hdu.write_key(fptr, "CDELT2", self.integration_time_s as f64 * 1000.0)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT2", "s")?;

            // SSINS-specific metadata
            hdu.write_key(fptr, "OBJECT", "EAVILS")?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;
            hdu.write_key(fptr, "POL", pol_names[pol_idx])?;

            // Add description
            hdu.write_key(
                fptr,
                "COMMENT",
                "EAVILS (Expected Amplitude of VisibILities Spectra)",
            )?;
        }

        // Write flag array
        let flag_dim = self.flag_array.dim();
        let flag_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[flag_dim.0, flag_dim.1],
        };
        let flag_extname = "EAVILS_FLAGS";
        let hdu = fptr.create_image(flag_extname, &flag_image_description)?;
        hdu.write_image(
            fptr,
            &self
                .flag_array
                .iter()
                .cloned()
                .map(|b| if b { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        )?;

        // Basic image info
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;

        // Time axis info
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz as f64)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz as f64)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // Frequency axis info
        hdu.write_key(fptr, "CTYPE2", "TIME")?;
        hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s as f64)?;
        hdu.write_key(fptr, "CDELT2", self.integration_time_s as f64 * 1000.0)?;
        hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT2", "s")?;

        // SSINS-specific metadata
        hdu.write_key(fptr, "OBJECT", "EAVILS")?;
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        // Add description
        hdu.write_key(
            fptr,
            "COMMENT",
            "EAVILS (Expected Amplitude of VisibILities Spectra)",
        )?;
        Ok(())
    }
}
