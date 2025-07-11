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
#[allow(clippy::upper_case_acronyms)]
pub(crate) struct SSINS {
    pub zscore: Array3<f32>,           // (times-1, frequencies, polarizations)
    pub diff_mean_amp_fp: Array2<f32>, // (frequencies, polarizations)
    pub flag_array: Array2<bool>,      // (times-1, frequencies)
    pub num_baselines: usize,
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
        // product is a 3D zscore array of shape (num_timesteps - 1, num_freqs, num_pols=4)
        // C = (4/pi - 1)
        // zscore = (N_bl / C).sqrt() * (mean_amp_tfp - mean_amp_fp) / mean_amp_fp

        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let mut diff_mean_amp_tfp = Array3::<f32>::zeros((num_timesteps - 1, num_freqs, 4));
        let mut diff_mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));

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

        diff_mean_amp_tfp /= num_baselines as f32;
        for t in 0..num_timesteps - 1 {
            for f in 0..num_freqs {
                for p in 0..4 {
                    diff_mean_amp_fp[[f, p]] += diff_mean_amp_tfp[[t, f, p]];
                }
            }
        }
        diff_mean_amp_fp /= (num_timesteps - 1) as f32;

        // subtract mean_amp_fp from mean_amp_tfp, divide by mean_amp_fp and multiply by mean_stdev_ratio
        let mut zscore = diff_mean_amp_tfp;
        let mean_stdev_ratio = (num_baselines as f32 / (4.0 / std::f32::consts::PI - 1.0)).sqrt();
        for t in 0..num_timesteps - 1 {
            for f in 0..num_freqs {
                for p in 0..4 {
                    zscore[[t, f, p]] = mean_stdev_ratio
                        * (zscore[[t, f, p]] - diff_mean_amp_fp[[f, p]])
                        / diff_mean_amp_fp[[f, p]];
                }
            }
        }
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
            zscore,
            diff_mean_amp_fp,
            flag_array,
            num_baselines,
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
        let imgset = amps_tfp_to_imageset(&aoflagger, self.zscore.view());
        let strategy_filename = strategy_filename
            .unwrap_or(aoflagger.FindStrategyFileGeneric(&String::from("minimal")));
        let flag_strategy = aoflagger.LoadStrategyFile(&strategy_filename.to_string());

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
        let pol_names = ["XX", "XY", "YX", "YY"];
        let (num_times, num_freqs, num_pols) = self.zscore.dim();

        for (pol_idx, &pol_name) in pol_names.iter().enumerate().take(num_pols) {
            let dim = [num_times, num_freqs];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("SSINS_POL={pol_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;

            // Write the data for this polarization
            let pol_data: Vec<f32> = self
                .zscore
                .slice(s![.., .., pol_idx])
                .iter()
                .copied()
                .collect();

            hdu.write_image(fptr, &pol_data)?;

            // Basic image info
            hdu.write_key(fptr, "BUNIT", "arbitrary")?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;

            // Time axis info
            hdu.write_key(fptr, "CTYPE1", "FREQ")?;
            hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
            hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "Hz")?;

            // Frequency axis info
            hdu.write_key(fptr, "CTYPE2", "TIME")?;
            hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
            hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT2", "s")?;

            // SSINS-specific metadata
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_BL", self.num_baselines as u32)?;
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
                .copied()
                .map(|b| if b { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        )?;

        // Basic image info
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;

        // Time axis info
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // Frequency axis info
        hdu.write_key(fptr, "CTYPE2", "TIME")?;
        hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
        hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
        hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT2", "s")?;

        // SSINS-specific metadata
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "SSINS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        // write diff_mean_amp_fp
        let diff_mean_amp_fp_dim = self.diff_mean_amp_fp.dim();
        let diff_mean_amp_fp_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[diff_mean_amp_fp_dim.0, diff_mean_amp_fp_dim.1],
        };
        let diff_mean_amp_fp_extname = "SSINS_DIFF_MEAN_AMP_FP";
        let hdu = fptr.create_image(
            diff_mean_amp_fp_extname,
            &diff_mean_amp_fp_image_description,
        )?;
        hdu.write_image(
            fptr,
            &self.diff_mean_amp_fp.iter().copied().collect::<Vec<_>>(),
        )?;
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;

        Ok(())
    }
}

#[allow(clippy::upper_case_acronyms)]
pub(crate) struct EAVILS {
    // Ntimes,Nbls,Nfreqs,Npols = uvd_or_data.shape
    // blmean_data = np.mean(np.abs(uvd_or_data),axis=1)
    // blmean_data_sub = blmean_data - np.mean(blmean_data,axis=0)
    // stdv_array = np.sqrt(np.mean(np.var(np.abs(uvd_or_data),axis=0,ddof=1),axis=0))
    // stdv_array =  stdv_array * np.full((Ntimes,Nfreqs,Npols),1) #giving it same shape as other arrays
    pub zscore: Array3<f32>, // (mean_amp_tfp - mean_amp_fp)/sqrt_mean_var_amp_fp (times, frequencies, polarizations)
    pub mean_amp_fp: Array2<f32>, // np.mean(blmean_data,axis=0) (frequencies, polarizations)
    // pub mean_amp_fbp: Array3<f32>, // np.mean(np.abs(uvd_or_data),axis=0) (frequencies, baselines, polarizations)
    // pub var_amp_fbp: Array3<f32>, // np.var(np.abs(uvd_or_data),axis=0,ddof=1) (frequencies, baselines, polarizations)
    pub sqrt_mean_var_amp_fp: Array2<f32>, // np.sqrt(np.mean(np.var(np.abs(uvd_or_data),axis=0,ddof=1),axis=0)) (frequencies, polarizations)
    pub flag_array: Array2<bool>,          // (times, frequencies)
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
        // mean_amp_fbp = jones_nodiff.sum(TIME) / num_timesteps
        // var_amp_fbp = (jones_nodiff.norm() - mean_amp_fbp).pow(2)
        // sqrt_mean_var_amp_fp = var_amp_fbp.mean(BL).sqrt()
        // zscore = (mean_amp_tfp - mean_amp_fp) / sqrt_mean_var_amp_fp

        let mut mean_amp_tfp = Array3::<f32>::zeros((num_timesteps, num_freqs, 4));
        let mut mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        let mut mean_amp_fbp = Array3::<f32>::zeros((num_freqs, num_baselines, 4));

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    let jones_nodiff = jones_array_tfb[[t, f, b]];
                    mean_amp_tfp[[t, f, 0]] += jones_nodiff[0].norm();
                    mean_amp_tfp[[t, f, 1]] += jones_nodiff[1].norm();
                    mean_amp_tfp[[t, f, 2]] += jones_nodiff[2].norm();
                    mean_amp_tfp[[t, f, 3]] += jones_nodiff[3].norm();
                    mean_amp_fbp[[f, b, 0]] += jones_nodiff[0].norm();
                    mean_amp_fbp[[f, b, 1]] += jones_nodiff[1].norm();
                    mean_amp_fbp[[f, b, 2]] += jones_nodiff[2].norm();
                    mean_amp_fbp[[f, b, 3]] += jones_nodiff[3].norm();
                }
            }
        }

        mean_amp_tfp /= num_baselines as f32;
        mean_amp_fbp /= num_timesteps as f32;

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for p in 0..4 {
                    mean_amp_fp[[f, p]] += mean_amp_tfp[[t, f, p]];
                }
            }
        }
        mean_amp_fp /= num_timesteps as f32;

        // calculate the time-variance of amplitude for each freq, bl, pol
        let mut var_amp_fbp = Array3::<f32>::zeros((num_freqs, num_baselines, 4));
        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    let jones_nodiff = jones_array_tfb[[t, f, b]];
                    var_amp_fbp[[f, b, 0]] +=
                        (jones_nodiff[0].norm() - mean_amp_fbp[[f, b, 0]]).powi(2);
                    var_amp_fbp[[f, b, 1]] +=
                        (jones_nodiff[1].norm() - mean_amp_fbp[[f, b, 1]]).powi(2);
                    var_amp_fbp[[f, b, 2]] +=
                        (jones_nodiff[2].norm() - mean_amp_fbp[[f, b, 2]]).powi(2);
                    var_amp_fbp[[f, b, 3]] +=
                        (jones_nodiff[3].norm() - mean_amp_fbp[[f, b, 3]]).powi(2);
                }
            }
        }

        var_amp_fbp /= num_timesteps as f32;

        // calculate mean_var_amp_fp
        let mut sqrt_mean_var_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        for f in 0..num_freqs {
            for b in 0..num_baselines {
                for p in 0..4 {
                    sqrt_mean_var_amp_fp[[f, p]] += var_amp_fbp[[f, b, p]];
                }
            }
        }
        sqrt_mean_var_amp_fp /= num_baselines as f32;
        sqrt_mean_var_amp_fp = sqrt_mean_var_amp_fp.sqrt();

        // calculate zscore
        let mut zscore = mean_amp_tfp;
        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for p in 0..4 {
                    zscore[[t, f, p]] =
                        (zscore[[t, f, p]] - mean_amp_fp[[f, p]]) / sqrt_mean_var_amp_fp[[f, p]];
                }
            }
        }

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
            zscore,
            mean_amp_fp,
            // mean_amp_fbp,
            // var_amp_fbp,
            sqrt_mean_var_amp_fp,
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
        let pol_names = ["XX", "XY", "YX", "YY"];
        let (num_times, num_freqs, num_pols) = self.zscore.dim();

        for (pol_idx, &pol_name) in pol_names.iter().enumerate().take(num_pols) {
            let dim = [num_times, num_freqs];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("EAVILS_POL={pol_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;

            // Write the data for this polarization
            let pol_data: Vec<f32> = self
                .zscore
                .slice(s![.., .., pol_idx])
                .iter()
                .copied()
                .collect();

            hdu.write_image(fptr, &pol_data)?;

            // Basic image info
            hdu.write_key(fptr, "BUNIT", "arbitrary")?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;

            // Time axis info
            hdu.write_key(fptr, "CTYPE1", "FREQ")?;
            hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
            hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "Hz")?;

            // Frequency axis info
            hdu.write_key(fptr, "CTYPE2", "TIME")?;
            hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
            hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT2", "s")?;

            // SSINS-specific metadata
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;
            hdu.write_key(fptr, "POL", pol_name)?;

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
                .copied()
                .map(|b| if b { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        )?;

        // Basic image info
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;

        // Time axis info
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // Frequency axis info
        hdu.write_key(fptr, "CTYPE2", "TIME")?;
        hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
        hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
        hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT2", "s")?;

        // SSINS-specific metadata
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        // Add description
        hdu.write_key(
            fptr,
            "COMMENT",
            "EAVILS (Expected Amplitude of VisibILities Spectra)",
        )?;

        // write mean_amp_fp
        let mean_amp_fp_dim = self.mean_amp_fp.dim();
        let mean_amp_fp_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[mean_amp_fp_dim.0, mean_amp_fp_dim.1],
        };
        let mean_amp_fp_extname = "EAVILS_MEAN_AMP_FP";
        let hdu = fptr.create_image(mean_amp_fp_extname, &mean_amp_fp_image_description)?;
        hdu.write_image(fptr, &self.mean_amp_fp.iter().copied().collect::<Vec<_>>())?;
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // write sqrt_mean_var_amp_fp
        let sqrt_mean_var_amp_fp_dim = self.sqrt_mean_var_amp_fp.dim();
        let sqrt_mean_var_amp_fp_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[sqrt_mean_var_amp_fp_dim.0, sqrt_mean_var_amp_fp_dim.1],
        };
        let sqrt_mean_var_amp_fp_extname = "EAVILS_SQRT_MEAN_VAR_AMP_FP";
        let hdu = fptr.create_image(
            sqrt_mean_var_amp_fp_extname,
            &sqrt_mean_var_amp_fp_image_description,
        )?;
        hdu.write_image(
            fptr,
            &self
                .sqrt_mean_var_amp_fp
                .iter()
                .copied()
                .collect::<Vec<_>>(),
        )?;
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        Ok(())
    }
}

pub(crate) struct AOFlagMetrics {
    pub occupancy_tf: Array2<f64>, // (times, frequencies)
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
}

impl AOFlagMetrics {
    pub(crate) fn new(
        flag_array_tfb: ArrayView3<bool>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = flag_array_tfb.dim();
        let mut occupancy_tf = Array2::<f64>::zeros((num_timesteps, num_freqs));

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    occupancy_tf[[t, f]] += flag_array_tfb[[t, f, b]] as u8 as f64;
                }
            }
        }

        let total_samples = num_timesteps * num_baselines;
        occupancy_tf /= total_samples as f64;

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
            occupancy_tf,
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
        let flag_dim: (usize, usize) = self.occupancy_tf.dim();
        let flag_image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[flag_dim.0, flag_dim.1],
        };
        let flag_extname = "AO_FLAG_METRICS";
        let hdu = fptr.create_image(flag_extname, &flag_image_description)?;
        hdu.write_image(fptr, &self.occupancy_tf.iter().copied().collect::<Vec<_>>())?;
        // Basic image info
        hdu.write_key(fptr, "BUNIT", "arbitrary")?;
        hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        hdu.write_key(fptr, "BZERO", 0.0f64)?;

        // Time axis info
        hdu.write_key(fptr, "CTYPE1", "FREQ")?;
        hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
        hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
        hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT1", "Hz")?;

        // Frequency axis info
        hdu.write_key(fptr, "CTYPE2", "TIME")?;
        hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
        hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
        hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
        hdu.write_key(fptr, "CUNIT2", "s")?;
        Ok(())
    }
}
