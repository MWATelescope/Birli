//! SSINS (Sky-Subtract Incoherent Noise Spectra)
//!
//! This module implements the SSINS algorithm, which is a method for
//! subtracting incoherent noise from the visibilities.
//!
//! The algorithm is described in:
//! - <https://arxiv.org/abs/1906.01093>
//! - <https://ssins.readthedocs.io/en/latest/>

use std::collections::HashMap;

use marlu::XyzGeocentric;

use crate::{
    delay_transform::{calculate_delay_channels, delay_transform, DelayTransformConfig},
    marlu::{
        fitsio::{
            images::{ImageDescription, ImageType},
            FitsFile,
        },
        mwalib::CorrelatorContext,
        ndarray::{s, Array2, Array3, Array4, ArrayView3},
        Jones, VisSelection,
    },
    FlagContext,
};

/// when you want to convert hyperdrive stokes order to standard stokes order.
/// fits standard:
/// −5 'XX' X parallel linear
/// −6 'YY' Y parallel linear
/// −7 'XY' XY cross linear
/// −8 'YX' YX cross linear
/// but hyperdrive is XX XY YX YY
pub fn hyperdrive_to_fits_stokes(pol: usize) -> usize {
    match pol {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 1,
        _ => panic!("Invalid pol: {}", pol),
    }
}

// Autocorrelation metrics
pub(crate) struct AutoMetrics {
    pub auto_sub_aptf: Array4<f32>, // auto with mean(time) subtracted (antenna, times, frequencies, polarizations)
    pub auto_spectrum_afp: Array3<f32>, // auto mean(time) (antenna, frequencies, polarizations)
    pub auto_delay_afp: Array3<f32>, // delay transform of auto mean(time) (antenna, delays, polarizations)
    pub antenna_names: Vec<String>,
    pub antenna_positions: Vec<XyzGeocentric>,
    pub antenna_ids: Vec<u32>,
    pub antenna_nums: Vec<u32>,
    pub rx_numbers: Vec<u32>,
    pub rx_slots: Vec<u32>,
    pub rx_types: Vec<String>,
    pub cable_flavours: Vec<String>,
    pub whitening_filters: Vec<bool>,

    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
}

impl AutoMetrics {
    pub(crate) fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());
        let sel_ant_pairs = chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        let sel_auto_pairs: HashMap<_, _> = sel_ant_pairs
            .iter()
            .enumerate()
            .filter(|&(_, (a, b))| a == b)
            .map(|(i, &(a, _))| (a, i))
            .collect();
        let num_sel_ants = sel_ant_pairs.iter().filter(|(a, b)| a == b).count();
        assert_eq!(num_sel_ants, sel_auto_pairs.len());
        let sel_auto_idxs = sel_auto_pairs.keys().copied().collect::<Vec<_>>();
        
        // Create a mapping from antenna numbers to their positions in the selected antennas list
        let mut ant_to_pos: HashMap<usize, usize> = HashMap::new();
        for (pos, &ant_num) in sel_auto_idxs.iter().enumerate() {
            ant_to_pos.insert(ant_num, pos);
        }
        
        let num_freqs = jones_array_tfb.dim().1;
        let num_timesteps = jones_array_tfb.dim().0;
        let num_pols = 4;
        let mut auto_power_aptf = Array4::<f32>::zeros((num_sel_ants, 4, num_timesteps, num_freqs));
        let mut auto_spectrum_afp = Array3::<f32>::zeros((num_sel_ants, num_freqs, 4));
        let num_unflagged_times_f32: f32 = timestep_flags.iter().filter(|&t| !t).count() as f32;
        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for (&a, &i) in &sel_auto_pairs {
                    let a_pos = ant_to_pos[&a]; // Get the position of antenna a in the selected antennas list
                    for p in 0..4 {
                        if chan_flags[f] {
                            auto_power_aptf[[a_pos, hyperdrive_to_fits_stokes(p), t, f]] = f32::NAN;
                            auto_spectrum_afp[[a_pos, f, hyperdrive_to_fits_stokes(p)]] = f32::NAN;
                        } else if timestep_flags[t] {
                            auto_power_aptf[[a_pos, hyperdrive_to_fits_stokes(p), t, f]] = f32::NAN;
                        } else {
                            let jones_p_norm = jones_array_tfb[[t, f, i]][p].norm();
                            auto_spectrum_afp[[a_pos, f, hyperdrive_to_fits_stokes(p)]] +=
                                jones_p_norm / num_unflagged_times_f32;
                            auto_power_aptf[[a_pos, hyperdrive_to_fits_stokes(p), t, f]] +=
                                jones_p_norm;
                        }
                    }
                }
            }
        }

        // take auto_power_aptf and subtract its mean(time), auto_spectrum_afp
        let mut auto_sub_aptf = auto_power_aptf;
        for a in 0..num_sel_ants {
            for p in 0..num_pols {
                for t in 0..num_timesteps {
                    for f in 0..num_freqs {
                        auto_sub_aptf[[a, hyperdrive_to_fits_stokes(p), t, f]] -=
                            auto_spectrum_afp[[a, f, hyperdrive_to_fits_stokes(p)]];
                    }
                }
            }
        }

        // variance over frequency for each antenna, time, and polarization
        // let mut auto_var_atp = Array3::<f32>::zeros((num_sel_ants, num_timesteps, num_pols));
        // for a in 0..num_sel_ants {
        //     for t in 0..num_timesteps {
        //         for p in 0..num_pols {
        //             auto_var_atp[[a, t, p]] = auto_sub_aptf
        //                 .slice(s![a, hyperdrive_to_fits_stokes(p), t, ..])
        //                 .var(0.0);
        //         }
        //     }
        // }

        let mut antenna_ids = Vec::<u32>::with_capacity(num_sel_ants);
        let mut antenna_names = Vec::<String>::with_capacity(num_sel_ants);
        let mut antenna_positions = Vec::<XyzGeocentric>::with_capacity(num_sel_ants);
        let mut antenna_nums = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_numbers = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_slots = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_types = Vec::<String>::with_capacity(num_sel_ants);
        let mut cable_flavours = Vec::<String>::with_capacity(num_sel_ants);
        let mut whitening_filters = Vec::<bool>::with_capacity(num_sel_ants);

        println!("sel_auto_idxs: {:?}", sel_auto_idxs);

        let freqs = corr_ctx.get_fine_chan_freqs_hz_array(
            &chunk_vis_sel.coarse_chan_range.clone().collect::<Vec<_>>(),
        );
        // use delay_transform to get the delay spectrum of each antenna in auto_spectrum_afp
        let delay_transform_config = DelayTransformConfig {
            min_delay_ns: 100.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1.0,
        };
        let freqs_arr = marlu::ndarray::Array1::from(freqs.clone());
        let delay_info = calculate_delay_channels(num_freqs, &freqs_arr, &delay_transform_config);
        let num_delays = delay_info.n_delay_channels;
        let mut auto_delay_afp = Array3::<f32>::zeros((num_sel_ants, num_delays, num_pols));
        for pol in 0..num_pols {
            // Convert spectrum to f64 and shape (ants, freqs)
            let spectrum_f64 = auto_spectrum_afp
                .slice(s![.., .., pol])
                .mapv(|x| x as f64)
                .to_owned();

            // Convert freqs to ndarray of f64
            let freqs_arr = marlu::ndarray::Array1::from(freqs.clone());
            let delay_spectrum =
                delay_transform(&spectrum_f64, &freqs_arr, &delay_transform_config)
                    .expect("delay_transform failed");
            auto_delay_afp
                .slice_mut(s![.., .., pol])
                .assign(&delay_spectrum.delay_spectrum.mapv(|x| x as f32));
        }

        corr_ctx
            .metafits_context
            .antennas
            .iter()
            .enumerate()
            .filter(|(i, _)| sel_auto_idxs.contains(i))
            .for_each(|(_, a)| {
                antenna_ids.push(a.ant);
                antenna_names.push(a.tile_name.clone());
                antenna_positions.push(XyzGeocentric {
                    x: a.north_m,
                    y: a.east_m,
                    z: a.height_m,
                });
                antenna_nums.push(a.tile_id);
                let rfinput_x = a.rfinput_x.clone();
                rx_numbers.push(rfinput_x.rec_number);
                rx_slots.push(rfinput_x.rec_slot_number);
                rx_types.push(rfinput_x.rec_type.to_string());
                cable_flavours.push(rfinput_x.flavour.clone());
                whitening_filters.push(rfinput_x.has_whitening_filter);
            });

        let timesteps = &corr_ctx.timesteps[chunk_vis_sel.timestep_range.clone()]
            .iter()
            .map(|ts| ts.gps_time_ms as f64 / 1000.0)
            .collect::<Vec<_>>();

        Self {
            auto_sub_aptf,
            auto_spectrum_afp,
            auto_delay_afp,
            antenna_names,
            antenna_positions,
            antenna_ids,
            antenna_nums,
            rx_numbers,
            rx_slots,
            rx_types,
            cable_flavours,
            whitening_filters,
            start_freq_hz: freqs[0],
            channel_width_hz: freqs[1] - freqs[0],
            start_time_gps_s: timesteps[0],
            integration_time_s: timesteps[1] - timesteps[0],
        }
    }

    pub(crate) fn save_to_fits(
        &self,
        fptr: &mut FitsFile,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let (num_ants, num_pols, _num_times, num_freqs) = self.auto_sub_aptf.dim();

        for a in 0..self.auto_sub_aptf.dim().0 {
            let antenna_name = self.antenna_names[a].clone();
            let dim = [
                self.auto_sub_aptf.dim().1,
                self.auto_sub_aptf.dim().2,
                self.auto_sub_aptf.dim().3,
            ];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("AUTO_SUB_ANT={antenna_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;
            hdu.write_image(
                fptr,
                &self
                    .auto_sub_aptf
                    .slice(s![a, .., .., ..])
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
            )?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;
            hdu.write_key(fptr, "CTYPE3", "FREQ")?;
            hdu.write_key(fptr, "CRVAL3", self.start_freq_hz)?;
            // hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
            hdu.write_key(fptr, "CRPIX3", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT3", "Hz")?;
            hdu.write_key(fptr, "CTYPE2", "TIME")?;
            hdu.write_key(fptr, "CRVAL2", self.start_time_gps_s)?;
            hdu.write_key(fptr, "CDELT2", self.integration_time_s)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;

            hdu.write_key(fptr, "CTYPE1", "STOKES")?;
            hdu.write_key(fptr, "CRVAL1", -5.0f64)?; // xx
            hdu.write_key(fptr, "CDELT1", -1.0f64)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;

            hdu.write_key(fptr, "ANTNAME", antenna_name)?;
            hdu.write_key(fptr, "ANT_ID", self.antenna_ids[a])?;
            hdu.write_key(fptr, "ANT_NUM", self.antenna_nums[a])?;
            hdu.write_key(fptr, "ANT_TYPE", self.rx_types[a].clone())?;
            hdu.write_key(fptr, "CABLE_FLAVOUR", self.cable_flavours[a].clone())?;
            hdu.write_key(fptr, "WHITENING_FILTER", self.whitening_filters[a] as i32)?;
            hdu.write_key(fptr, "RX_NUMBER", self.rx_numbers[a])?;
            hdu.write_key(fptr, "RX_SLOT", self.rx_slots[a])?;
            hdu.write_key(fptr, "RX_TYPE", self.rx_types[a].clone())?;

            let position = self.antenna_positions[a];
            hdu.write_key(fptr, "OBSGEO-X", position.x)?;
            hdu.write_key(fptr, "OBSGEO-Y", position.y)?;
            hdu.write_key(fptr, "OBSGEO-Z", position.z)?;
        }

        for pol_idx in 0..num_pols {
            // this is in fits standard order, not hyperdrive order
            let pol_name = ["XX", "YY", "XY", "YX"][pol_idx];
            let dim = [num_ants, num_freqs];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("AUTO_POL={pol_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;
            hdu.write_image(
                fptr,
                &self
                    .auto_spectrum_afp
                    .slice(s![.., .., pol_idx])
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
            )?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;
            hdu.write_key(fptr, "CTYPE1", "FREQ")?;
            hdu.write_key(fptr, "CRVAL1", self.start_freq_hz)?;
            // if left uncommented, carta will not display it properly.
            // hdu.write_key(fptr, "CDELT1", self.channel_width_hz)?;
            hdu.write_key(fptr, "CHAN_WIDTH", self.channel_width_hz)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "Hz")?;
            hdu.write_key(fptr, "CTYPE2", "BASELINE")?;
            hdu.write_key(fptr, "CRVAL2", 0.0f64)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
        }

        // // write out auto_var_atp
        // for pol_idx in 0..num_pols {
        //     let pol_name = ["XX", "YY", "XY", "YX"][pol_idx];
        //     let dim = [num_ants, num_timesteps];
        //     let image_description = ImageDescription {
        //         data_type: ImageType::Double,
        //         dimensions: &dim,
        //     };
        //     let extname = format!("AUTO_VAR_POL={pol_name}");
        //     let hdu = fptr.create_image(&extname, &image_description)?;
        //     hdu.write_image(
        //         fptr,
        //         &self
        //             .auto_var_atp
        //             .slice(s![.., .., pol_idx])
        //             .iter()
        //             .copied()
        //             .collect::<Vec<_>>(),
        //     )?;
        //     hdu.write_key(fptr, "BSCALE", 1.0f64)?;
        //     hdu.write_key(fptr, "BZERO", 0.0f64)?;
        //     hdu.write_key(fptr, "CTYPE1", "TIME")?;
        //     hdu.write_key(fptr, "CRVAL1", self.start_time_gps_s)?;
        //     hdu.write_key(fptr, "CDELT1", self.integration_time_s)?;
        //     hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
        //     hdu.write_key(fptr, "CUNIT1", "s")?;
        //     hdu.write_key(fptr, "CTYPE2", "BASELINE")?;
        // }

        let num_delays = self.auto_delay_afp.dim().1;
        // write out auto_delay_spectrum_afp
        for pol_idx in 0..num_pols {
            let pol_name = ["XX", "YY", "XY", "YX"][pol_idx];
            let dim = [num_ants, num_delays];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("AUTO_DELAY_POL={pol_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;
            hdu.write_image(
                fptr,
                &self
                    .auto_delay_afp
                    .slice(s![.., .., pol_idx])
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
            )?;
            hdu.write_key(fptr, "BSCALE", 1.0f64)?;
            hdu.write_key(fptr, "BZERO", 0.0f64)?;
            hdu.write_key(fptr, "CTYPE1", "DELAY")?;
            hdu.write_key(fptr, "CRVAL1", 0.0f64)?;
            hdu.write_key(fptr, "CRPIX1", 1.0f64)?;
            hdu.write_key(fptr, "CUNIT1", "ns")?;
            hdu.write_key(fptr, "CTYPE2", "BASELINE")?;
            hdu.write_key(fptr, "CRVAL2", 0.0f64)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
        }

        Ok(())
    }
}

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
        flag_ctx: &FlagContext,
    ) -> Self {
        // incoherently averaged (along baseline) difference between the visibilities in time.
        // product is a 3D zscore array of shape (num_timesteps - 1, num_freqs, num_pols=4)
        // C = (4/pi - 1)
        // zscore = (N_bl / C).sqrt() * (mean_amp_tfp - mean_amp_fp) / mean_amp_fp

        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();

        if num_timesteps < 2 {
            panic!("SSINS requires at least 2 timesteps");
        }
        let mut diff_mean_amp_tfp = Array3::<f32>::zeros((num_timesteps - 1, num_freqs, 4));
        let mut diff_mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));

        let flag_array = Array2::<bool>::default((num_timesteps - 1, num_freqs));
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let num_unflagged_diff_timesteps = timestep_flags
            .iter()
            .zip(timestep_flags[1..].iter())
            .filter(|&(a, b)| !a && !b)
            .count();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        for t in 0..num_timesteps - 1 {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    for p in 0..4 {
                        if chan_flags[f] {
                            diff_mean_amp_tfp[[t, f, p]] = f32::NAN;
                            diff_mean_amp_fp[[f, p]] = f32::NAN;
                        } else if timestep_flags[t] || timestep_flags[t + 1] {
                            diff_mean_amp_tfp[[t, f, p]] = f32::NAN;
                        } else {
                            let jones_diff_p_norm = (jones_array_tfb[[t, f, b]]
                                - jones_array_tfb[[t + 1, f, b]])[p]
                                .norm();
                            diff_mean_amp_tfp[[t, f, p]] += jones_diff_p_norm;
                            diff_mean_amp_fp[[f, p]] += jones_diff_p_norm;
                        }
                    }
                }
            }
        }

        if num_baselines > 0 {
            diff_mean_amp_tfp /= num_baselines as f32;
            diff_mean_amp_fp /= num_baselines as f32;
        }
        if num_unflagged_diff_timesteps > 0 {
            diff_mean_amp_fp /= num_unflagged_diff_timesteps as f32;
        }

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
        let integration_time_s = timesteps_nodiff[1] - timesteps_nodiff[0];
        // Compute the average time between adjacent timesteps
        let timesteps_diff: Vec<f64> = timesteps_nodiff[1..]
            .iter()
            .zip(&timesteps_nodiff[..num_timesteps - 1])
            .map(|(a, b)| (a + b) / 2.0)
            .collect();

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
        flag_ctx: &FlagContext,
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        // mean_amp_fbp = jones_nodiff.sum(TIME) / num_timesteps
        // var_amp_fbp = (jones_nodiff.norm() - mean_amp_fbp).pow(2)
        // sqrt_mean_var_amp_fp = var_amp_fbp.mean(BL).sqrt()
        // zscore = (mean_amp_tfp - mean_amp_fp) / sqrt_mean_var_amp_fp

        let mut mean_amp_tfp = Array3::<f32>::zeros((num_timesteps, num_freqs, 4));
        // EAVILS_MEAN_AMP_FP
        let mut mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        let mut mean_amp_fbp = Array3::<f32>::zeros((num_freqs, num_baselines, 4));
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let num_unflagged_timesteps = timestep_flags.iter().filter(|&t| !t).count();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    for p in 0..4 {
                        if chan_flags[f] {
                            mean_amp_fbp[[f, b, p]] = f32::NAN;
                            mean_amp_tfp[[t, f, p]] = f32::NAN;
                        } else if timestep_flags[t] {
                            mean_amp_tfp[[t, f, p]] = f32::NAN;
                        } else {
                            let jones_p_norm = jones_array_tfb[[t, f, b]][p].norm();
                            mean_amp_tfp[[t, f, p]] += jones_p_norm;
                            mean_amp_fbp[[f, b, p]] += jones_p_norm;
                        }
                    }
                }
            }
        }

        mean_amp_tfp /= num_baselines as f32;
        if num_unflagged_timesteps > 0 {
            mean_amp_fbp /= num_unflagged_timesteps as f32;
        }

        for t in 0..num_timesteps {
            if timestep_flags[t] {
                continue;
            }
            for f in 0..num_freqs {
                for p in 0..4 {
                    if chan_flags[f] {
                        mean_amp_fp[[f, p]] = f32::NAN;
                        continue;
                    }
                    mean_amp_fp[[f, p]] += mean_amp_tfp[[t, f, p]];
                }
            }
        }
        mean_amp_fp /= num_unflagged_timesteps as f32;

        // calculate the time-variance of amplitude for each freq, bl, pol
        // variance is mean of the squared deviations from the mean
        let mut var_amp_fbp = Array3::<f32>::zeros((num_freqs, num_baselines, 4));
        for t in 0..num_timesteps {
            if timestep_flags[t] {
                continue;
            }
            for f in 0..num_freqs {
                for b in 0..num_baselines {
                    for p in 0..4 {
                        if chan_flags[f] {
                            var_amp_fbp[[f, b, p]] = f32::NAN;
                            continue;
                        }
                        let jones_nodiff = jones_array_tfb[[t, f, b]];
                        var_amp_fbp[[f, b, p]] +=
                            (jones_nodiff[p].norm() - mean_amp_fbp[[f, b, p]]).powi(2);
                    }
                }
            }
        }

        if num_unflagged_timesteps > 0 {
            var_amp_fbp /= num_unflagged_timesteps as f32;
        }

        // calculate mean_var_amp_fp, EAVILS_SQRT_MEAN_VAR_AMP_FP
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

        // calculate zscore (EAVILS_POL)
        let mut zscore = mean_amp_tfp;
        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for p in 0..4 {
                    zscore[[t, f, p]] =
                        (zscore[[t, f, p]] - mean_amp_fp[[f, p]]) / sqrt_mean_var_amp_fp[[f, p]];
                }
            }
        }
        // multiply zscore by sqrt_num_unflagged_baselines
        zscore *= (num_baselines as f32).sqrt();

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

#[cfg(feature = "aoflagger")]
pub(crate) struct AOFlagMetrics {
    pub occupancy_tf: Array2<f64>, // (times, frequencies)
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
}

#[cfg(feature = "aoflagger")]
impl AOFlagMetrics {
    /// Compute occupancy across baselines excluding autocorrelations
    ///
    /// occupancy[t,f] = fraction of cross-correlation baselines flagged at (t,f)
    fn compute_cross_occupancy(
        flag_array_tfb: ArrayView3<bool>,
        timestep_flags: &[bool],
        chan_flags: &[bool],
        ant_pairs: &[(usize, usize)],
    ) -> Array2<f64> {
        let (num_timesteps, num_freqs, _num_baselines) = flag_array_tfb.dim();
        let mut occupancy_tf = Array2::<f64>::zeros((num_timesteps, num_freqs));

        // indices of cross-correlation baselines (exclude a==b)
        let cross_baseline_indices: Vec<usize> = ant_pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &(a, b))| if a != b { Some(i) } else { None })
            .collect();
        let num_cross_baselines = cross_baseline_indices.len();

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                if timestep_flags[t] || chan_flags[f] {
                    occupancy_tf[[t, f]] = f64::NAN;
                    continue;
                }
                let mut sum_flags = 0u32;
                for &b_idx in &cross_baseline_indices {
                    sum_flags += flag_array_tfb[[t, f, b_idx]] as u32;
                }
                if num_cross_baselines > 0 {
                    occupancy_tf[[t, f]] = sum_flags as f64 / num_cross_baselines as f64;
                } else {
                    occupancy_tf[[t, f]] = f64::NAN;
                }
            }
        }

        occupancy_tf
    }
    pub(crate) fn new(
        flag_array_tfb: ArrayView3<bool>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());
        let sel_ant_pairs = chunk_vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        let occupancy_tf = Self::compute_cross_occupancy(
            flag_array_tfb,
            &timestep_flags,
            &chan_flags,
            &sel_ant_pairs,
        );

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

#[cfg(all(test, feature = "aoflagger"))]
mod aoflagmetrics_tests {
    use super::AOFlagMetrics;
    use crate::marlu::ndarray::{Array2, Array3};

    #[test]
    fn test_crosses_only_occupancy_helper() {
        // Two timesteps, two freqs, three baselines: (0,0) auto, (0,1) cross, (1,1) auto
        let mut flags = Array3::<bool>::from_elem((2, 2, 3), false);
        // Set autos flagged everywhere
        for t in 0..2 {
            for f in 0..2 {
                flags[[t, f, 0]] = true; // auto 0-0
                flags[[t, f, 2]] = true; // auto 1-1
            }
        }
        // Cross baseline initially unflagged
        flags[[0, 0, 1]] = false;
        flags[[0, 1, 1]] = false;
        flags[[1, 0, 1]] = true; // flag cross at (t=1,f=0)
        flags[[1, 1, 1]] = false;

        let timestep_flags = vec![false, false];
        let chan_flags = vec![false, false];
        let ant_pairs = vec![(0usize, 0usize), (0usize, 1usize), (1usize, 1usize)];

        let occ = AOFlagMetrics::compute_cross_occupancy(
            flags.view(),
            &timestep_flags,
            &chan_flags,
            &ant_pairs,
        );

        // Expected: fraction over crosses only. Only one cross baseline exists.
        let mut expected = Array2::<f64>::zeros((2, 2));
        expected[[0, 0]] = 0.0; // cross false
        expected[[0, 1]] = 0.0; // cross false
        expected[[1, 0]] = 1.0; // cross true
        expected[[1, 1]] = 0.0; // cross false

        for t in 0..2 {
            for f in 0..2 {
                assert!(
                    (occ[[t, f]] - expected[[t, f]]).abs() < 1e-12,
                    "mismatch at (t={}, f={}): got {}, expected {}",
                    t,
                    f,
                    occ[[t, f]],
                    expected[[t, f]]
                );
            }
        }
    }
}

#[cfg(test)]
mod autometrics_tests {
    use super::AutoMetrics;
    use crate::{
        marlu::{
            mwalib::CorrelatorContext,
            ndarray::Array3,
            Jones, VisSelection,
        },
        FlagContext,
    };
    use tempfile::tempdir;

    #[test]
    fn test_autometrics_with_non_sequential_antenna_selection() {
        // This test specifically verifies the fix for the array indexing bug
        // where antenna numbers (like 1, 2) were used directly as array indices
        // instead of their positions in the selected antennas list (0, 1)
        
        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        // Create correlator context
        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        
        // Create visibility selection with non-sequential antenna selection (1, 2)
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[1, 2]);
        
        // Verify that we have the expected antenna pairs
        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        assert_eq!(ant_pairs, vec![(1, 1), (1, 2), (2, 2)]);
        
        // Create flag context
        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        
        // Create a small jones array for testing (2 timesteps, 1 frequency, 3 baselines)
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 3));
        
        // This should not panic with the fixed indexing
        let auto_metrics = AutoMetrics::new(
            jones_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
        );
        
        // Verify that the metrics were created successfully
        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 2); // 2 selected antennas
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 2); // 2 selected antennas
        assert_eq!(auto_metrics.auto_delay_afp.dim().0, 2); // 2 selected antennas
        
        // Verify antenna names and IDs are correct
        assert_eq!(auto_metrics.antenna_names.len(), 2);
        assert_eq!(auto_metrics.antenna_ids.len(), 2);
        assert_eq!(auto_metrics.antenna_nums.len(), 2);
    }

    #[test]
    fn test_autometrics_with_sequential_antenna_selection() {
        // Test with sequential antenna selection (0, 1) to ensure
        // the fix doesn't break the normal case
        
        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[0, 1]);
        
        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        assert_eq!(ant_pairs, vec![(0, 0), (0, 1), (1, 1)]);
        
        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 3));
        
        let auto_metrics = AutoMetrics::new(
            jones_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
        );
        
        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 2);
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 2);
        assert_eq!(auto_metrics.auto_delay_afp.dim().0, 2);
    }

    #[test]
    fn test_autometrics_with_single_antenna_selection() {
        // Test with single antenna selection to ensure edge case works
        
        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[5]);
        
        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        assert_eq!(ant_pairs, vec![(5, 5)]);
        
        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 1));
        
        let auto_metrics = AutoMetrics::new(
            jones_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
        );
        
        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 1);
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 1);
        assert_eq!(auto_metrics.auto_delay_afp.dim().0, 1);
    }

    #[test]
    fn test_autometrics_save_to_fits() {
        // Test that AutoMetrics can be saved to FITS without errors
        
        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths = vec![
            "tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits",
        ];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
        
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[1, 2]);
        
        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 3));
        
        let auto_metrics = AutoMetrics::new(
            jones_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
        );
        
        // Test saving to FITS
        let tmp_dir = tempdir().unwrap();
        let fits_path = tmp_dir.path().join("test_autometrics.fits");
        
        let mut fptr = crate::marlu::fitsio::FitsFile::create(&fits_path).open().unwrap();
        assert!(auto_metrics.save_to_fits(&mut fptr).is_ok());
        
        // Verify file was created
        assert!(fits_path.exists());
        assert!(fits_path.metadata().unwrap().len() > 0);
    }
}
