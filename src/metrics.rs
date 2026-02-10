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
        ndarray::{s, Array1, Array2, Array3, Array4, ArrayView3, Axis},
        num_complex::Complex,
        Jones, VisSelection,
    },
    FlagContext,
};

/// Minimal metadata required for metrics calculation
/// This allows metrics to work with any data source, not just mwalib CorrelatorContext
#[derive(Debug, Clone)]
pub struct MetricsContext {
    /// Antenna information for selected antennas
    pub antennas: Vec<AntennaMetadata>,
    /// Fine channel frequencies in Hz
    pub fine_chan_freqs_hz: Vec<f64>,
    /// GPS timestamps in seconds
    pub timestamps_s: Vec<f64>,
    /// Antenna pairs (indices into the antennas vec) - corresponds to baseline dimension
    pub antenna_pairs: Vec<(usize, usize)>,
}

/// Antenna metadata needed for metrics
#[derive(Debug, Clone)]
pub struct AntennaMetadata {
    pub tile_name: String,
    pub tile_id: u32,
    pub ant_id: u32,
    /// ENU position in metres
    pub east_m: f64,
    pub north_m: f64,
    pub height_m: f64,
    /// Receiver information
    pub rec_number: u32,
    pub rec_slot_number: u32,
    pub rec_type: String,
    pub cable_flavour: String,
    pub has_whitening_filter: bool,
}

impl MetricsContext {
    /// Create MetricsContext from mwalib CorrelatorContext and VisSelection
    /// This is a convenience function for backward compatibility
    pub fn from_mwalib(corr_ctx: &CorrelatorContext, vis_sel: &VisSelection) -> Self {
        let fine_chan_freqs_hz = corr_ctx
            .get_fine_chan_freqs_hz_array(&vis_sel.coarse_chan_range.clone().collect::<Vec<_>>());

        let timestamps_s: Vec<f64> = corr_ctx.timesteps[vis_sel.timestep_range.clone()]
            .iter()
            .map(|ts| ts.gps_time_ms as f64 / 1000.0)
            .collect();

        let antenna_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);

        // Get unique antennas from antenna_pairs
        let mut unique_ant_indices = std::collections::HashSet::new();
        for &(a, b) in &antenna_pairs {
            unique_ant_indices.insert(a);
            unique_ant_indices.insert(b);
        }
        let mut sorted_ant_indices: Vec<usize> = unique_ant_indices.into_iter().collect();
        sorted_ant_indices.sort_unstable();

        let antennas: Vec<AntennaMetadata> = sorted_ant_indices
            .iter()
            .map(|&idx| {
                let ant = &corr_ctx.metafits_context.antennas[idx];
                AntennaMetadata {
                    tile_name: ant.tile_name.clone(),
                    tile_id: ant.tile_id,
                    ant_id: ant.ant,
                    east_m: ant.east_m,
                    north_m: ant.north_m,
                    height_m: ant.height_m,
                    rec_number: ant.rfinput_x.rec_number,
                    rec_slot_number: ant.rfinput_x.rec_slot_number,
                    rec_type: ant.rfinput_x.rec_type.to_string(),
                    cable_flavour: ant.rfinput_x.flavour.clone(),
                    has_whitening_filter: ant.rfinput_x.has_whitening_filter,
                }
            })
            .collect();

        Self {
            antennas,
            fine_chan_freqs_hz,
            timestamps_s,
            antenna_pairs,
        }
    }
}

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
pub struct AutoMetrics {
    pub auto_sub_aptf: Array4<f32>, // auto with mean(time) subtracted (antenna, times, frequencies, polarizations)
    pub auto_spectrum_afp: Array3<f32>, // auto mean(time) (antenna, frequencies, polarizations)
    pub auto_coeffs_apo: Array3<f32>, // polynomial coeffs (antenna, pol, order) - raw frequency basis
    pub auto_delay_afp: Array3<f32>, // delay transform of auto mean(time) (antenna, delays, polarizations)
    pub antenna_names: Vec<String>,
    pub antenna_ids: Vec<u32>,
    pub antenna_nums: Vec<u32>,
    pub rx_numbers: Vec<u32>,
    pub rx_slots: Vec<u32>,
    pub rx_types: Vec<String>,
    pub cable_flavours: Vec<String>,
    pub whitening_filters: Vec<bool>,
    pub antenna_positions: Vec<XyzGeocentric>,

    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
}

impl AutoMetrics {
    // RUST_LOG=birli=debug cargo run --release -- --sel-ants 5 4 20 19 --provided-chan-ranges --flag-init 0 --metrics-out metrics_1119683928.fits -m tests/data/1119683928_picket/1119683928.metafits tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits 2>&1 | tee birli.log

    /// Legacy constructor using CorrelatorContext (for backward compatibility)
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let metadata = MetricsContext::from_mwalib(corr_ctx, chunk_vis_sel);
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        Self::new_from_metadata(jones_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create AutoMetrics from visibility data and metadata
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();

        // Map from baseline index to auto-correlation antenna index in metadata.antennas
        let mut ant_idx_to_metadata_pos: HashMap<usize, usize> = HashMap::new();
        for (pos, ant_meta) in metadata.antennas.iter().enumerate() {
            // Find this antenna in antenna_pairs to get its original index
            for &(a_idx, b_idx) in &metadata.antenna_pairs {
                if a_idx == b_idx {
                    // This is an auto, check if it matches
                    // We need to map back - assume antenna_pairs uses indices that match position in original antenna list
                    // For now, we'll use a different approach
                    ant_idx_to_metadata_pos.insert(pos, pos);
                    break;
                }
            }
        }

        // Find auto-correlations in the baseline list
        let sel_auto_pairs: HashMap<usize, usize> = metadata
            .antenna_pairs
            .iter()
            .enumerate()
            .filter(|&(_, &(a, b))| a == b)
            .map(|(bl_idx, &(a, _))| (a, bl_idx))
            .collect();

        let num_sel_ants = sel_auto_pairs.len();
        let sel_auto_idxs = {
            let mut sel_auto_idxs = sel_auto_pairs.keys().copied().collect::<Vec<_>>();
            sel_auto_idxs.sort_unstable();
            sel_auto_idxs
        };

        // Create a mapping from antenna indices to their positions in the selected antennas list
        let mut ant_to_pos: HashMap<usize, usize> = HashMap::new();
        for (pos, &ant_num) in sel_auto_idxs.iter().enumerate() {
            ant_to_pos.insert(ant_num, pos);
        }

        let num_pols = 4;
        let mut auto_power_aptf = Array4::<f32>::zeros((num_sel_ants, 4, num_timesteps, num_freqs));
        let mut auto_spectrum_afp = Array3::<f32>::zeros((num_sel_ants, num_freqs, 4));
        let num_unflagged_times_f32: f32 = timestep_flags.iter().filter(|&t| !t).count() as f32;
        if (num_unflagged_times_f32 - 0.0).abs() < 0.00001 {
            panic!("all timesteps are flagged");
        }
        if (chan_flags.iter().filter(|&c| !c).count()) == 0 {
            panic!("all channels are flagged");
        }

        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                for (&a, &i) in &sel_auto_pairs {
                    let a_pos = ant_to_pos[&a];
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
        // fit polynomials (in frequency, Hz) to auto_spectrum_afp for each antenna and polarization
        let poly_order = 2;
        let freqs = &metadata.fine_chan_freqs_hz;
        let mut auto_coeffs_apo = Array3::<f32>::zeros((num_sel_ants, num_pols, poly_order + 1));

        for a in 0..num_sel_ants {
            for p in 0..num_pols {
                let spectrum = auto_spectrum_afp.slice(s![a, .., p]);
                let y: Vec<f32> = spectrum.iter().copied().collect();

                if let Some(coeffs) = crate::math::fit_polynomial(&freqs, &y, poly_order) {
                    for (k, &c) in coeffs.iter().enumerate() {
                        auto_coeffs_apo[[a, p, k]] = c;
                    }
                } else {
                    auto_coeffs_apo.slice_mut(s![a, p, ..]).fill(f32::NAN);
                }
            }
        }

        // Subtract mean
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

        let mut antenna_ids = Vec::<u32>::with_capacity(num_sel_ants);
        let mut antenna_names = Vec::<String>::with_capacity(num_sel_ants);
        let mut antenna_positions = Vec::<XyzGeocentric>::with_capacity(num_sel_ants);
        let mut antenna_nums = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_numbers = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_slots = Vec::<u32>::with_capacity(num_sel_ants);
        let mut rx_types = Vec::<String>::with_capacity(num_sel_ants);
        let mut cable_flavours = Vec::<String>::with_capacity(num_sel_ants);
        let mut whitening_filters = Vec::<bool>::with_capacity(num_sel_ants);

        // Delay transform
        let delay_transform_config = DelayTransformConfig {
            min_delay_ns: 100.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1.0,
        };
        let freqs_arr = marlu::ndarray::Array1::from(metadata.fine_chan_freqs_hz.clone());
        let delay_info = calculate_delay_channels(num_freqs, &freqs_arr, &delay_transform_config);
        let num_delays = delay_info.n_delay_channels;
        let mut auto_delay_afp = Array3::<f32>::zeros((num_sel_ants, num_delays, num_pols));

        for pol in 0..num_pols {
            let spectrum_f64 = auto_spectrum_afp
                .slice(s![.., .., pol])
                .mapv(|x| x as f64)
                .to_owned();

            let delay_spectrum =
                delay_transform(&spectrum_f64, &freqs_arr, &delay_transform_config)
                    .expect("delay_transform failed");
            auto_delay_afp
                .slice_mut(s![.., .., pol])
                .assign(&delay_spectrum.delay_spectrum.mapv(|x| x as f32));
        }

        // Collect antenna metadata for selected autos (sel_auto_idxs are metafits indices; metadata.antennas is 0..n by same order)
        for (pos, _) in sel_auto_idxs.iter().enumerate() {
            let ant = &metadata.antennas[pos];
            antenna_ids.push(ant.ant_id);
            antenna_names.push(ant.tile_name.clone());
            antenna_positions.push(XyzGeocentric {
                x: ant.north_m,
                y: ant.east_m,
                z: ant.height_m,
            });
            antenna_nums.push(ant.tile_id);
            rx_numbers.push(ant.rec_number);
            rx_slots.push(ant.rec_slot_number);
            rx_types.push(ant.rec_type.clone());
            cable_flavours.push(ant.cable_flavour.clone());
            whitening_filters.push(ant.has_whitening_filter);
        }

        Self {
            auto_sub_aptf,
            auto_spectrum_afp,
            auto_coeffs_apo,
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
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: if metadata.fine_chan_freqs_hz.len() > 1 {
                metadata.fine_chan_freqs_hz[1] - metadata.fine_chan_freqs_hz[0]
            } else {
                0.0
            },
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: if metadata.timestamps_s.len() > 1 {
                metadata.timestamps_s[1] - metadata.timestamps_s[0]
            } else {
                1.0 // Default to 1 second if only one timestep
            },
        }
    }

    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> Result<(), Box<dyn std::error::Error>> {
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

        // Write polynomial coefficients to binary table
        let extname = "AUTO_COEFFS";
        let poly_order = self.auto_coeffs_apo.dim().2 - 1;
        let num_coeffs = poly_order + 1;

        let mut ant_names_col = Vec::with_capacity(num_ants);
        let mut ant_ids_col = Vec::with_capacity(num_ants);
        let mut ant_nums_col = Vec::with_capacity(num_ants);
        let mut ant_types_col = Vec::with_capacity(num_ants);
        let mut cable_flavours_col = Vec::with_capacity(num_ants);
        let mut whitening_filters_col = Vec::with_capacity(num_ants);
        let mut rx_numbers_col = Vec::with_capacity(num_ants);
        let mut rx_slots_col = Vec::with_capacity(num_ants);
        let mut rx_types_col = Vec::with_capacity(num_ants);
        let mut obsgeo_x_col = Vec::with_capacity(num_ants);
        let mut obsgeo_y_col = Vec::with_capacity(num_ants);
        let mut obsgeo_z_col = Vec::with_capacity(num_ants);

        let mut coeffs_xx_col = Vec::with_capacity(num_ants * num_coeffs);
        let mut coeffs_yy_col = Vec::with_capacity(num_ants * num_coeffs);
        let mut coeffs_xy_col = Vec::with_capacity(num_ants * num_coeffs);
        let mut coeffs_yx_col = Vec::with_capacity(num_ants * num_coeffs);

        // hyperdrive order: XX XY YX YY
        // fits order: XX YY XY YX (we want to match what we write in other extensions)
        // But auto_coeffs_apo is indexed by hyperdrive_to_fits_stokes(p) if we were consistent?
        // Wait, auto_spectrum_afp is indexed by hyperdrive_to_fits_stokes(p).
        // Let's check AutoMetrics::new.
        //   auto_spectrum_afp[[a_pos, f, hyperdrive_to_fits_stokes(p)]]
        // So axis 2 of auto_spectrum_afp is fits stokes order.
        // Therefore axis 1 of auto_coeffs_apo is fits stokes order.
        // FITS Order: 0:XX, 1:YY, 2:XY, 3:YX

        for a in 0..num_ants {
            ant_names_col.push(self.antenna_names[a].clone());
            ant_ids_col.push(self.antenna_ids[a]);
            ant_nums_col.push(self.antenna_nums[a]);
            ant_types_col.push(self.rx_types[a].clone());
            cable_flavours_col.push(self.cable_flavours[a].clone());
            whitening_filters_col.push(self.whitening_filters[a] as i32);
            rx_numbers_col.push(self.rx_numbers[a]);
            rx_slots_col.push(self.rx_slots[a]);
            rx_types_col.push(self.rx_types[a].clone());

            let position = self.antenna_positions[a];
            obsgeo_x_col.push(position.x);
            obsgeo_y_col.push(position.y);
            obsgeo_z_col.push(position.z);

            // XX (index 0)
            for k in 0..num_coeffs {
                coeffs_xx_col.push(self.auto_coeffs_apo[[a, 0, k]] as f64);
            }
            // YY (index 1)
            for k in 0..num_coeffs {
                coeffs_yy_col.push(self.auto_coeffs_apo[[a, 1, k]] as f64);
            }
            // XY (index 2)
            for k in 0..num_coeffs {
                coeffs_xy_col.push(self.auto_coeffs_apo[[a, 2, k]] as f64);
            }
            // YX (index 3)
            for k in 0..num_coeffs {
                coeffs_yx_col.push(self.auto_coeffs_apo[[a, 3, k]] as f64);
            }
        }

        use crate::marlu::fitsio::tables::{ColumnDataType, ColumnDescription};

        let columns = vec![
            ColumnDescription::new("ANT_NAME")
                .with_type(ColumnDataType::String)
                .that_repeats(32) // String width
                .create()?,
            ColumnDescription::new("ANT_ID")
                .with_type(ColumnDataType::Int)
                .create()?,
            ColumnDescription::new("ANT_NUM")
                .with_type(ColumnDataType::Int)
                .create()?,
            ColumnDescription::new("ANT_TYPE")
                .with_type(ColumnDataType::String)
                .that_repeats(32)
                .create()?,
            ColumnDescription::new("CABLE_FLAVOUR")
                .with_type(ColumnDataType::String)
                .that_repeats(32)
                .create()?,
            ColumnDescription::new("WHITENING_FILTER")
                .with_type(ColumnDataType::Int)
                .create()?,
            ColumnDescription::new("RX_NUMBER")
                .with_type(ColumnDataType::Int)
                .create()?,
            ColumnDescription::new("RX_SLOT")
                .with_type(ColumnDataType::Int)
                .create()?,
            ColumnDescription::new("RX_TYPE")
                .with_type(ColumnDataType::String)
                .that_repeats(32)
                .create()?,
            ColumnDescription::new("OBSGEO-X")
                .with_type(ColumnDataType::Double)
                .create()?,
            ColumnDescription::new("OBSGEO-Y")
                .with_type(ColumnDataType::Double)
                .create()?,
            ColumnDescription::new("OBSGEO-Z")
                .with_type(ColumnDataType::Double)
                .create()?,
            ColumnDescription::new("COEFFS_XX")
                .with_type(ColumnDataType::Double)
                .that_repeats(num_coeffs)
                .create()?,
            ColumnDescription::new("COEFFS_YY")
                .with_type(ColumnDataType::Double)
                .that_repeats(num_coeffs)
                .create()?,
            ColumnDescription::new("COEFFS_XY")
                .with_type(ColumnDataType::Double)
                .that_repeats(num_coeffs)
                .create()?,
            ColumnDescription::new("COEFFS_YX")
                .with_type(ColumnDataType::Double)
                .that_repeats(num_coeffs)
                .create()?,
        ];

        let mut table_hdu = fptr.create_table(extname, &columns)?;

        table_hdu.write_col(fptr, "ANT_NAME", &ant_names_col)?;
        table_hdu.write_col(fptr, "ANT_ID", &ant_ids_col)?;
        table_hdu.write_col(fptr, "ANT_NUM", &ant_nums_col)?;
        table_hdu.write_col(fptr, "ANT_TYPE", &ant_types_col)?;
        table_hdu.write_col(fptr, "CABLE_FLAVOUR", &cable_flavours_col)?;
        table_hdu.write_col(fptr, "WHITENING_FILTER", &whitening_filters_col)?;
        table_hdu.write_col(fptr, "RX_NUMBER", &rx_numbers_col)?;
        table_hdu.write_col(fptr, "RX_SLOT", &rx_slots_col)?;
        table_hdu.write_col(fptr, "RX_TYPE", &rx_types_col)?;
        table_hdu.write_col(fptr, "OBSGEO-X", &obsgeo_x_col)?;
        table_hdu.write_col(fptr, "OBSGEO-Y", &obsgeo_y_col)?;
        table_hdu.write_col(fptr, "OBSGEO-Z", &obsgeo_z_col)?;

        table_hdu.write_col(fptr, "COEFFS_XX", &coeffs_xx_col)?;
        table_hdu.write_col(fptr, "COEFFS_YY", &coeffs_yy_col)?;
        table_hdu.write_col(fptr, "COEFFS_XY", &coeffs_xy_col)?;
        table_hdu.write_col(fptr, "COEFFS_YX", &coeffs_yx_col)?;

        table_hdu.write_key(fptr, "TELESCOP", "MWA")?;
        table_hdu.write_key(fptr, "POLY_ORD", poly_order as u32)?;

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

/// SSINS (Sky-Subtract Incoherent Noise Spectra)
///
/// incoherently averaged (along baseline) difference between the visibilities in time.
/// product is a 3D zscore array of shape (num_timesteps - 1, num_freqs, num_pols=4)
/// C = (4/pi - 1)
/// zscore = (N_bl / C).sqrt() * (mean_amp_tfp - mean_amp_fp) / mean_amp_fp
#[allow(clippy::upper_case_acronyms)]
pub struct SSINS {
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
    /// Legacy constructor using CorrelatorContext (for backward compatibility)
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let metadata = MetricsContext::from_mwalib(corr_ctx, chunk_vis_sel);
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        Self::new_from_metadata(jones_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create SSINS from visibility data and metadata
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();

        if num_timesteps < 2 {
            panic!("SSINS requires at least 2 timesteps");
        }
        let mut diff_mean_amp_tfp = Array3::<f32>::zeros((num_timesteps - 1, num_freqs, 4));
        let mut diff_mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));

        let flag_array = Array2::<bool>::default((num_timesteps - 1, num_freqs));
        let num_unflagged_diff_timesteps = timestep_flags
            .iter()
            .zip(timestep_flags[1..].iter())
            .filter(|&(a, b)| !a && !b)
            .count();

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

        let integration_time_s = if metadata.timestamps_s.len() > 1 {
            metadata.timestamps_s[1] - metadata.timestamps_s[0]
        } else {
            1.0 // Default to 1 second if only one timestep
        };
        // Compute the average time between adjacent timesteps
        let timesteps_diff: Vec<f64> = if num_timesteps > 1 {
            metadata.timestamps_s[1..]
                .iter()
                .zip(&metadata.timestamps_s[..num_timesteps - 1])
                .map(|(a, b)| (a + b) / 2.0)
                .collect()
        } else {
            vec![metadata.timestamps_s[0]]
        };

        let freq_width_hz = metadata.fine_chan_freqs_hz[1] - metadata.fine_chan_freqs_hz[0];

        Self {
            zscore,
            diff_mean_amp_fp,
            flag_array,
            num_baselines,
            start_time_gps_s: timesteps_diff[0],
            integration_time_s,
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: freq_width_hz,
        }
    }

    #[cfg(feature = "aoflagger")]
    pub fn flag(&mut self, strategy_filename: Option<String>) {
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

    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> Result<(), Box<dyn std::error::Error>> {
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
pub struct EAVILS {
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
    /// Legacy constructor using CorrelatorContext (for backward compatibility)
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let metadata = MetricsContext::from_mwalib(corr_ctx, chunk_vis_sel);
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        Self::new_from_metadata(jones_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create EAVILS from visibility data and metadata
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();

        let mut mean_amp_tfp = Array3::<f32>::zeros((num_timesteps, num_freqs, 4));
        let mut mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        let mut mean_amp_fbp = Array3::<f32>::zeros((num_freqs, num_baselines, 4));
        let num_unflagged_timesteps = timestep_flags.iter().filter(|&t| !t).count();

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
        zscore *= (num_baselines as f32).sqrt();

        let flag_array = Array2::<bool>::default((num_timesteps, num_freqs));

        let integration_time_s = if metadata.timestamps_s.len() > 1 {
            metadata.timestamps_s[1] - metadata.timestamps_s[0]
        } else {
            1.0 // Default to 1 second if only one timestep
        };
        let freq_width_hz = if metadata.fine_chan_freqs_hz.len() > 1 {
            metadata.fine_chan_freqs_hz[1] - metadata.fine_chan_freqs_hz[0]
        } else {
            1.0 // Default to 1 Hz if only one frequency
        };

        Self {
            zscore,
            mean_amp_fp,
            sqrt_mean_var_amp_fp,
            flag_array,
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s,
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: freq_width_hz,
        }
    }

    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> Result<(), Box<dyn std::error::Error>> {
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
pub struct AOFlagMetrics {
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

    /// Legacy constructor using CorrelatorContext (for backward compatibility)
    pub fn new(
        flag_array_tfb: ArrayView3<bool>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let metadata = MetricsContext::from_mwalib(corr_ctx, chunk_vis_sel);
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        Self::new_from_metadata(flag_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create AOFlagMetrics from flag data and metadata
    pub fn new_from_metadata(
        flag_array_tfb: ArrayView3<bool>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let occupancy_tf = Self::compute_cross_occupancy(
            flag_array_tfb,
            timestep_flags,
            chan_flags,
            &metadata.antenna_pairs,
        );

        let integration_time_s = if metadata.timestamps_s.len() > 1 {
            metadata.timestamps_s[1] - metadata.timestamps_s[0]
        } else {
            1.0 // Default to 1 second if only one timestep
        };
        let freq_width_hz = if metadata.fine_chan_freqs_hz.len() > 1 {
            metadata.fine_chan_freqs_hz[1] - metadata.fine_chan_freqs_hz[0]
        } else {
            1.0 // Default to 1 Hz if only one frequency
        };

        Self {
            occupancy_tf,
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s,
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: freq_width_hz,
        }
    }

    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> Result<(), Box<dyn std::error::Error>> {
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

#[cfg(feature = "aoflagger")]
pub struct CrossMetrics {
    pub short_ant_delay_pol_adp: Array3<f32>, // (antennas, delay bin, polarizations)
    pub antenna_names: Vec<String>,
    pub antenna_ids: Vec<u32>,
    pub antenna_nums: Vec<u32>,
    pub start_freq_hz: f64,
    pub channel_width_hz: f64,
    pub start_time_gps_s: f64,
    pub integration_time_s: f64,
}

#[cfg(feature = "aoflagger")]
impl CrossMetrics {
    /// Legacy constructor using CorrelatorContext (for backward compatibility)
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
        baseline_cutoff_m: f32,
    ) -> Self {
        let metadata = MetricsContext::from_mwalib(corr_ctx, chunk_vis_sel);
        let timestep_flags = flag_ctx.timestep_flags[chunk_vis_sel.timestep_range.clone()].to_vec();
        let chan_flags = flag_ctx.get_raw_chan_flags(&chunk_vis_sel.coarse_chan_range.clone());

        Self::new_from_metadata(
            jones_array_tfb,
            &metadata,
            &timestep_flags,
            &chan_flags,
            baseline_cutoff_m,
        )
    }

    /// Create CrossMetrics from visibility data and metadata
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
        baseline_cutoff_m: f32,
    ) -> Self {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();

        // Identify selected antennas
        let mut sel_ants_set = std::collections::HashSet::new();
        for &(a, b) in &metadata.antenna_pairs {
            sel_ants_set.insert(a);
            sel_ants_set.insert(b);
        }
        let mut sel_ants_sorted: Vec<usize> = sel_ants_set.into_iter().collect();
        sel_ants_sorted.sort_unstable();

        let num_sel_ants = sel_ants_sorted.len();
        let mut ant_to_pos = HashMap::new();
        for (i, &ant_idx) in sel_ants_sorted.iter().enumerate() {
            ant_to_pos.insert(ant_idx, i);
        }

        let mut antenna_names = Vec::with_capacity(num_sel_ants);
        let mut antenna_ids = Vec::with_capacity(num_sel_ants);
        let mut antenna_nums = Vec::with_capacity(num_sel_ants);

        for (pos, _) in sel_ants_sorted.iter().enumerate() {
            let ant = &metadata.antennas[pos];
            antenna_names.push(ant.tile_name.clone());
            antenna_ids.push(ant.ant_id);
            antenna_nums.push(ant.tile_id);
        }

        // Setup delay transform
        let freqs = &metadata.fine_chan_freqs_hz;
        let delay_transform_config = DelayTransformConfig {
            min_delay_ns: 100.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1.0,
        };
        let freqs_arr = marlu::ndarray::Array1::from(freqs.clone());
        let delay_info = calculate_delay_channels(num_freqs, &freqs_arr, &delay_transform_config);
        let num_delays = delay_info.n_delay_channels;

        let mut short_ant_delay_pol_adp = Array3::<f32>::zeros((num_sel_ants, num_delays, 4));
        let mut ant_baseline_counts = Array2::<u32>::zeros((num_sel_ants, 4));

        // Iterate over baselines
        for (b_idx, &(a_idx, b_idx_ant)) in metadata.antenna_pairs.iter().enumerate() {
            if a_idx == b_idx_ant {
                continue; // Skip autos
            }

            let pos_a = ant_to_pos[&a_idx];
            let pos_b = ant_to_pos[&b_idx_ant];
            let ant_a = &metadata.antennas[pos_a];
            let ant_b = &metadata.antennas[pos_b];

            let dx = ant_a.north_m - ant_b.north_m;
            let dy = ant_a.east_m - ant_b.east_m;
            let dz = ant_a.height_m - ant_b.height_m;
            let len_sq = dx * dx + dy * dy + dz * dz;

            if len_sq > (baseline_cutoff_m * baseline_cutoff_m) as f64 {
                continue;
            }

            // Calculate complex spectrum for this baseline
            // Mean over time of complex V
            let mut complex_spectrum = Array2::<Complex<f64>>::zeros((num_freqs, 4));
            let mut counts = Array2::<u32>::zeros((num_freqs, 4));

            for t in 0..num_timesteps {
                if timestep_flags[t] {
                    continue;
                }
                for f in 0..num_freqs {
                    if chan_flags[f] {
                        continue;
                    }

                    for p in 0..4 {
                        let val = jones_array_tfb[[t, f, b_idx]][p];
                        if val.re.is_finite() && val.im.is_finite() {
                            complex_spectrum[[f, p]] += Complex::new(val.re as f64, val.im as f64);
                            counts[[f, p]] += 1;
                        }
                    }
                }
            }

            // Normalize complex_spectrum
            for f in 0..num_freqs {
                for p in 0..4 {
                    if counts[[f, p]] > 0 {
                        complex_spectrum[[f, p]] /= counts[[f, p]] as f64;
                    }
                }
            }

            // Perform delay transform for each pol
            for p in 0..4 {
                // Create spectrum array for this pol (magnitude for delay transform)
                let spec_col: Array1<f64> = complex_spectrum
                    .column(p)
                    .iter()
                    .map(|c| c.norm())
                    .collect();

                // If we have valid data
                if spec_col.iter().any(|x| *x > 0.0) {
                    if let Ok(delay_res) = delay_transform(
                        &spec_col.insert_axis(Axis(0)),
                        &freqs_arr,
                        &delay_transform_config,
                    ) {
                        let delay_spec = delay_res.delay_spectrum.row(0);

                        // Accumulate to antennas
                        if let Some(&idx_a) = ant_to_pos.get(&a_idx) {
                            let mut slice = short_ant_delay_pol_adp.slice_mut(s![idx_a, .., p]);
                            for d in 0..num_delays {
                                slice[d] += delay_spec[d] as f32;
                            }
                            ant_baseline_counts[[idx_a, p]] += 1;
                        }
                        if let Some(&idx_b) = ant_to_pos.get(&b_idx_ant) {
                            let mut slice = short_ant_delay_pol_adp.slice_mut(s![idx_b, .., p]);
                            for d in 0..num_delays {
                                slice[d] += delay_spec[d] as f32;
                            }
                            ant_baseline_counts[[idx_b, p]] += 1;
                        }
                    }
                }
            }
        }

        // Normalize by number of baselines per antenna
        for a in 0..num_sel_ants {
            for p in 0..4 {
                let count = ant_baseline_counts[[a, p]];
                if count > 0 {
                    let mut slice = short_ant_delay_pol_adp.slice_mut(s![a, .., p]);
                    for d in 0..num_delays {
                        slice[d] /= count as f32;
                    }
                } else {
                    let mut slice = short_ant_delay_pol_adp.slice_mut(s![a, .., p]);
                    slice.fill(f32::NAN);
                }
            }
        }

        Self {
            short_ant_delay_pol_adp,
            antenna_names,
            antenna_ids,
            antenna_nums,
            start_freq_hz: freqs[0],
            channel_width_hz: freqs[1] - freqs[0],
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: if metadata.timestamps_s.len() > 1 {
                metadata.timestamps_s[1] - metadata.timestamps_s[0]
            } else {
                0.0
            },
        }
    }

    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> Result<(), Box<dyn std::error::Error>> {
        let (num_ants, num_delays, num_pols) = self.short_ant_delay_pol_adp.dim();

        for pol_idx in 0..num_pols {
            let pol_name = ["XX", "YY", "XY", "YX"][pol_idx];
            let dim = [num_ants, num_delays];
            let image_description = ImageDescription {
                data_type: ImageType::Double,
                dimensions: &dim,
            };
            let extname = format!("CROSS_DELAY_POL={pol_name}");
            let hdu = fptr.create_image(&extname, &image_description)?;

            hdu.write_image(
                fptr,
                &self
                    .short_ant_delay_pol_adp
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
            hdu.write_key(fptr, "CTYPE2", "ANTENNA")?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "CROSS_METRICS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;

            // Write antenna names/IDs if possible as table or keywords?
            // AutoMetrics writes separate HDU per antenna for "AUTO_SUB_ANT", but "AUTO_POL" is (ants, freqs).
            // Here we have (ants, delays).
            // It is simpler to write one image per pol.
        }
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
        marlu::{mwalib::CorrelatorContext, ndarray::Array3, Jones, VisSelection},
        FlagContext,
    };
    use tempfile::tempdir;

    #[test]
    fn test_autometrics_with_non_sequential_antenna_selection() {
        // This test specifically verifies the fix for the array indexing bug
        // where antenna numbers (like 1, 2) were used directly as array indices
        // instead of their positions in the selected antennas list (0, 1)

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths =
            vec!["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"];

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
        let auto_metrics = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        // Verify that the metrics were created successfully
        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 2); // 2 selected antennas
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 2); // 2 selected antennas
                                                               // auto_coeffs_apo is (ants, pols, coeffs).
                                                               // 2 ants, 4 pols, 4 coeffs (order 3)
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().0, 2);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().1, 4);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().2, 3);
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
        let gpufits_paths =
            vec!["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[0, 1]);

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        assert_eq!(ant_pairs, vec![(0, 0), (0, 1), (1, 1)]);

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 3));

        let auto_metrics = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 2);
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 2);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().0, 2);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().1, 4);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().2, 3);
        assert_eq!(auto_metrics.auto_delay_afp.dim().0, 2);
    }

    #[test]
    fn test_autometrics_with_single_antenna_selection() {
        // Test with single antenna selection to ensure edge case works

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths =
            vec!["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[5]);

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        assert_eq!(ant_pairs, vec![(5, 5)]);

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 1));

        let auto_metrics = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        assert_eq!(auto_metrics.auto_sub_aptf.dim().0, 1);
        assert_eq!(auto_metrics.auto_spectrum_afp.dim().0, 1);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().0, 1);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().1, 4);
        assert_eq!(auto_metrics.auto_coeffs_apo.dim().2, 3);
        assert_eq!(auto_metrics.auto_delay_afp.dim().0, 1);
    }

    #[test]
    fn test_autometrics_save_to_fits() {
        // Test that AutoMetrics can be saved to FITS without errors

        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths =
            vec!["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[1, 2]);

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let jones_array = Array3::<Jones<f32>>::zeros((2, 1, 3));

        let auto_metrics = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        // Test saving to FITS
        let tmp_dir = tempdir().unwrap();
        let fits_path = tmp_dir.path().join("test_autometrics.fits");

        let mut fptr = crate::marlu::fitsio::FitsFile::create(&fits_path)
            .open()
            .unwrap();
        assert!(auto_metrics.save_to_fits(&mut fptr).is_ok());

        // Verify file was created
        assert!(fits_path.exists());
        assert!(fits_path.metadata().unwrap().len() > 0);
    }
}
#[cfg(all(test, feature = "aoflagger"))]
mod crossmetrics_tests {
    use super::CrossMetrics;
    use crate::{
        marlu::{mwalib::CorrelatorContext, ndarray::Array3, Jones, VisSelection},
        FlagContext,
    };
    #[test]
    fn test_crossmetrics_new() {
        let metafits_path = "tests/data/1119683928_picket/1119683928.metafits";
        let gpufits_paths =
            vec!["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"];

        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        // Use sequential antenna selection for simplicity (0, 1)
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[0, 1]);

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        // Should have (0,0), (0,1), (1,1) if both are in input

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        // Create dummy jones array (times, freqs, baselines)
        // 3 baselines: 0-0, 0-1, 1-1. We need to match vis_sel structure.
        let num_baselines = ant_pairs.len();
        let (num_timesteps, num_freqs) = (2, 4);
        let mut jones_array =
            Array3::<Jones<f32>>::zeros((num_timesteps, num_freqs, num_baselines));

        // Populate cross baseline (index 1) with some value
        // Assuming 0-1 is at index 1
        for t in 0..num_timesteps {
            for f in 0..num_freqs {
                jones_array[[t, f, 1]] = Jones::identity();
            }
        }

        let baseline_cutoff_m = 10000.0; // Large enough to include baselines

        let cross_metrics = CrossMetrics::new(
            jones_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
            baseline_cutoff_m,
        );

        assert_eq!(cross_metrics.short_ant_delay_pol_adp.dim().0, 2); // 2 antennas
                                                                      // Check values are not NaN (at least for valid pols)
                                                                      // Since we put constant value, delay transform should have peak at 0 delay (DC component).
    }
}
