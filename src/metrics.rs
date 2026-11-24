//! Quality-assurance metrics for preprocessed visibilities.
//!
//! These metrics are cheap summaries of visibilities, which are computed
//! after the corrections and flagging, before averaging. They are
//! written to a FITS file with `birli --metrics-out <PATH>`, or they can be
//! computed and written one at a time with this module.
//!
//! | metric | what it is | needs |
//! | --- | --- | --- |
//! | [`AutoMetrics`] | per-tile autocorrelation bandpass, its polynomial fit and delay transform, and the autocorrelations with the bandpass subtracted | autocorrelations |
//! | [`SSINS`] | z-scores of the incoherent noise spectrum: the baseline-averaged amplitude of the time-differenced cross-correlations, as in [SSINS](https://ssins.readthedocs.io/en/latest/) ([Wilensky et al. 2019](https://arxiv.org/abs/1906.01093)) | at least 2 timesteps |
//! | [`EAVILS`] | z-scores of the baseline-averaged amplitude of the visibilities, without time differencing | |
//! | `AOFlagMetrics` | fraction of cross-correlation baselines flagged at each time and frequency | `aoflagger` feature |
//! | [`CrossMetrics`] | per-tile delay transform of the cross-correlations on short baselines | |
//!
//! Every metric has three ways to make it:
//! - `new` takes the [`mwalib::CorrelatorContext`](CorrelatorContext),
//!   [`VisSelection`] and [`FlagContext`] of an observation that was read with
//!   Birli.
//! - `new_from_metadata` takes a [`MetricsContext`], so that the metrics can be
//!   computed for visibilities from any source.
//! - a builder (like [`AutoMetricsBuilder`]) takes the visibilities in chunks
//!   of timesteps, for observations that don't fit in memory. The metrics
//!   are the same however the timesteps are chunked.
//!
//! [`MetricsWriter`] drives all of the builders, and writes every metric to
//! one file. This is what `birli --metrics-out` uses.
//!
//! Every metric also has a `save_to_fits` method, which appends the metric's
//! HDUs to an open FITS file. The layout of these HDUs is described in
//! [`doc/metrics.md`](https://github.com/MWATelescope/Birli/blob/main/doc/metrics.md).
//!
//! # Conventions
//!
//! - Visibility and flag arrays have dimensions `[timestep][channel][baseline]`,
//!   like the arrays in the rest of Birli.
//! - Timesteps in `timestep_flags` and channels in `chan_flags` are left out of
//!   every metric. Only [`SSINS`] and `AOFlagMetrics` look at the flags of
//!   individual visibilities.
//! - [`Jones`] matrices are in the order `XX, XY, YX, YY`. [`SSINS`] and
//!   [`EAVILS`] keep this order. [`AutoMetrics`] arrays are in the order of the
//!   FITS Stokes axis, `XX, YY, XY, YX`, see [`hyperdrive_to_fits_stokes`].
//!
//! # Examples
//!
//! ```rust
//! use birli::{
//!     marlu::{fitsio::FitsFile, mwalib::CorrelatorContext},
//!     metrics::{AutoMetrics, EAVILS, SSINS},
//!     FlagContext, VisSelection,
//! };
//! use tempfile::tempdir;
//!
//! // a legacy observation with 4 timesteps of 2 coarse channels, each with 2
//! // fine channels
//! let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
//! let gpufits_paths = vec![
//!     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
//!     "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
//!     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
//!     "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
//! ];
//! let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
//!
//! // read the first 4 tiles
//! let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
//! vis_sel.retain_antennas(&corr_ctx.metafits_context, &[0, 1, 2, 3]);
//! let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
//! let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
//! let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
//! let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
//! flag_ctx
//!     .set_flags(
//!         flag_array.view_mut(),
//!         &vis_sel.timestep_range,
//!         &vis_sel.coarse_chan_range,
//!         &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
//!     )
//!     .unwrap();
//! birli::io::read_mwalib(&vis_sel, &corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
//!     .unwrap();
//!
//! let auto_metrics = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);
//! let ssins = SSINS::new(jones_array.view(), flag_array.view(), &corr_ctx, &vis_sel, &flag_ctx);
//! let eavils = EAVILS::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);
//!
//! // one bandpass per tile and polarisation
//! assert_eq!(auto_metrics.auto_spectrum_afp.dim(), (4, 4, 4));
//! // one fewer timestep than the data, because SSINS differences in time
//! assert_eq!(ssins.zscore.dim(), (3, 4, 4));
//! assert_eq!(eavils.zscore.dim(), (4, 4, 4));
//!
//! let tmp_dir = tempdir().unwrap();
//! let mut fptr = FitsFile::create(tmp_dir.path().join("metrics.fits")).open().unwrap();
//! auto_metrics.save_to_fits(&mut fptr).unwrap();
//! ssins.save_to_fits(&mut fptr).unwrap();
//! eavils.save_to_fits(&mut fptr).unwrap();
//! ```

mod writer;
pub use writer::MetricsWriter;

use std::collections::BTreeSet;

use rayon::prelude::*;

use marlu::ENH;

use crate::{
    delay_transform::{calculate_delay_channels, delay_transform, DelayTransformConfig},
    marlu::{
        fitsio::{
            hdu::FitsHdu,
            images::{ImageDescription, ImageType, WriteImage},
            tables::{ColumnDataType, ColumnDescription},
            FitsFile,
        },
        mwalib::CorrelatorContext,
        ndarray::{
            s, Array1, Array2, Array3, Array4, ArrayView2, ArrayView3, ArrayViewMut1,
            ArrayViewMut2, Axis,
        },
        num_complex::Complex,
        Jones, VisSelection,
    },
    math::fit_polynomial,
    FlagContext,
};

/// The result of writing to a FITS file.
type FitsResult<T> = Result<T, crate::marlu::fitsio::errors::Error>;

/// Names of the polarisations in the order of a [`Jones`] matrix.
const JONES_POL_NAMES: [&str; 4] = ["XX", "XY", "YX", "YY"];

/// Names of the polarisations in the order of the FITS Stokes axis (-5, -6,
/// -7, -8), see [`hyperdrive_to_fits_stokes`].
const FITS_POL_NAMES: [&str; 4] = ["XX", "YY", "XY", "YX"];

/// The delays in the delay transforms of [`AutoMetrics`] and [`CrossMetrics`].
const DELAY_TRANSFORM_CONFIG: DelayTransformConfig = DelayTransformConfig {
    min_delay_ns: 100.0,
    max_delay_ns: 3000.0,
    target_delay_res_ns: 1.0,
};

/// The order of the polynomial that [`AutoMetrics`] fits to each bandpass.
const AUTO_POLY_ORDER: usize = 2;

/// The metadata that is needed to calculate metrics.
///
/// This allows metrics to work with visibilities from any source, not just
/// those that were read with an mwalib [`CorrelatorContext`].
#[derive(Debug, Clone)]
pub struct MetricsContext {
    /// Information about each antenna that is in `antenna_pairs`, in order of
    /// increasing antenna index. So if `antenna_pairs` only has antennas 2 and
    /// 5, this has two elements: antenna 2, then antenna 5.
    pub antennas: Vec<AntennaMetadata>,
    /// Centre frequency of each fine channel \[Hz\]
    pub fine_chan_freqs_hz: Vec<f64>,
    /// GPS time of each timestep \[s\]
    pub timestamps_s: Vec<f64>,
    /// The antenna indices of each baseline, in the order of the baseline axis
    /// of the visibilities. Autocorrelations have the same index twice.
    pub antenna_pairs: Vec<(usize, usize)>,
}

/// Antenna metadata needed for metrics
#[derive(Debug, Clone)]
pub struct AntennaMetadata {
    /// Tile name from the metafits
    pub tile_name: String,
    /// Tile ID from the metafits
    pub tile_id: u32,
    /// Antenna index
    pub ant_id: u32,
    /// East position in metres (ENU)
    pub east_m: f64,
    /// North position in metres (ENU)
    pub north_m: f64,
    /// Height in metres (ENU)
    pub height_m: f64,
    /// Receiver number
    pub rec_number: u32,
    /// Slot number of this tile in its receiver
    pub rec_slot_number: u32,
    /// Receiver type
    pub rec_type: String,
    /// Cable flavour
    pub cable_flavour: String,
    /// Whether the cable has a whitening filter
    pub has_whitening_filter: bool,
}

impl MetricsContext {
    /// Create a [`MetricsContext`] for the timesteps, coarse channels and
    /// baselines of `vis_sel` in the observation `corr_ctx`.
    ///
    /// The receiver and cable of a tile are those of its X polarisation.
    ///
    /// # Panics
    ///
    /// Panics if `vis_sel` is not a selection of `corr_ctx`.
    pub fn from_mwalib(corr_ctx: &CorrelatorContext, vis_sel: &VisSelection) -> Self {
        let fine_chan_freqs_hz = corr_ctx
            .get_fine_chan_freqs_hz_array(&vis_sel.coarse_chan_range.clone().collect::<Vec<_>>())
            .expect("failed to get fine channel frequencies");

        let timestamps_s: Vec<f64> = corr_ctx.timesteps[vis_sel.timestep_range.clone()]
            .iter()
            .map(|ts| ts.gps_time_ms as f64 / 1000.0)
            .collect();

        let antenna_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);

        let antennas: Vec<AntennaMetadata> = selected_antennas(&antenna_pairs)
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

    /// The time between the first two timesteps \[s\], or 1 if there is only
    /// one timestep.
    fn integration_time_s(&self) -> f64 {
        if self.timestamps_s.len() > 1 {
            self.timestamps_s[1] - self.timestamps_s[0]
        } else {
            1.0
        }
    }

    /// The frequency difference between the first two channels \[Hz\], or 1 if
    /// there is only one channel.
    fn channel_width_hz(&self) -> f64 {
        if self.fine_chan_freqs_hz.len() > 1 {
            self.fine_chan_freqs_hz[1] - self.fine_chan_freqs_hz[0]
        } else {
            1.0
        }
    }
}

/// The sorted indices of the antennas which are in any of `antenna_pairs`.
/// [`MetricsContext::antennas`] is in this order.
fn selected_antennas(antenna_pairs: &[(usize, usize)]) -> Vec<usize> {
    antenna_pairs
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// The metadata, timestep flags and channel flags for the selection `vis_sel`
/// of the observation `corr_ctx`.
fn selection_metadata(
    corr_ctx: &CorrelatorContext,
    vis_sel: &VisSelection,
    flag_ctx: &FlagContext,
) -> (MetricsContext, Vec<bool>, Vec<bool>) {
    (
        MetricsContext::from_mwalib(corr_ctx, vis_sel),
        flag_ctx.timestep_flags[vis_sel.timestep_range.clone()].to_vec(),
        flag_ctx.get_raw_chan_flags(&vis_sel.coarse_chan_range),
    )
}

/// Convert a polarisation index in the order of a [`Jones`] matrix (as in
/// hyperdrive: `XX, XY, YX, YY`) to an index in the order of the FITS Stokes
/// axis:
///
/// | Stokes | polarisation | FITS index | Jones index |
/// | --- | --- | --- | --- |
/// | -5 | `XX` | 0 | 0 |
/// | -6 | `YY` | 1 | 3 |
/// | -7 | `XY` | 2 | 1 |
/// | -8 | `YX` | 3 | 2 |
///
/// # Panics
///
/// Panics if `pol` is not less than 4.
pub fn hyperdrive_to_fits_stokes(pol: usize) -> usize {
    match pol {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 1,
        _ => panic!("Invalid pol: {pol}"),
    }
}

/// Append an image HDU called `extname` to a FITS file.
///
/// `dimensions` are in the order of the array that `data` was flattened from,
/// which is the reverse of the order of the FITS axes: the last dimension is
/// FITS axis 1.
fn write_image_hdu<T: WriteImage + Copy>(
    fptr: &mut FitsFile,
    extname: &str,
    data_type: ImageType,
    dimensions: &[usize],
    data: impl IntoIterator<Item = T>,
) -> FitsResult<FitsHdu> {
    let hdu = fptr.create_image(
        extname,
        &ImageDescription {
            data_type,
            dimensions,
        },
    )?;
    hdu.write_image(fptr, &data.into_iter().collect::<Vec<_>>())?;
    hdu.write_key(fptr, "BSCALE", 1.0f64)?;
    hdu.write_key(fptr, "BZERO", 0.0f64)?;
    Ok(hdu)
}

/// A linear axis of a FITS image.
struct FitsAxis<'a> {
    /// What the axis is (`CTYPEn`)
    ctype: &'a str,
    /// The value at the first pixel (`CRVALn`)
    crval: f64,
    /// The increment between pixels (`CDELTn`), if it is written.
    cdelt: Option<f64>,
    /// The unit of the values (`CUNITn`), if it is written.
    cunit: Option<&'a str>,
}

impl FitsAxis<'_> {
    /// Describe FITS axis number `axis` (starting from 1) of `hdu`.
    fn write(&self, fptr: &mut FitsFile, hdu: &FitsHdu, axis: usize) -> FitsResult<()> {
        hdu.write_key(fptr, &format!("CTYPE{axis}"), self.ctype)?;
        hdu.write_key(fptr, &format!("CRVAL{axis}"), self.crval)?;
        if let Some(cdelt) = self.cdelt {
            hdu.write_key(fptr, &format!("CDELT{axis}"), cdelt)?;
        }
        hdu.write_key(fptr, &format!("CRPIX{axis}"), 1.0f64)?;
        if let Some(cunit) = self.cunit {
            hdu.write_key(fptr, &format!("CUNIT{axis}"), cunit)?;
        }
        Ok(())
    }
}

/// The time and frequency axes of a metric.
#[derive(Debug, Clone, Copy)]
struct TimeFreqAxes {
    /// GPS time of the first timestep \[s\]
    start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    integration_time_s: f64,
    /// Centre frequency of the first channel \[Hz\]
    start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    channel_width_hz: f64,
}

impl TimeFreqAxes {
    const fn freq(&self) -> FitsAxis<'static> {
        FitsAxis {
            ctype: "FREQ",
            crval: self.start_freq_hz,
            cdelt: Some(self.channel_width_hz),
            cunit: Some("Hz"),
        }
    }

    const fn time(&self) -> FitsAxis<'static> {
        FitsAxis {
            ctype: "TIME",
            crval: self.start_time_gps_s,
            cdelt: Some(self.integration_time_s),
            cunit: Some("s"),
        }
    }

    /// Describe an image with dimensions `[timestep][channel]`: frequency is
    /// FITS axis 1 and time is FITS axis 2.
    fn write(&self, fptr: &mut FitsFile, hdu: &FitsHdu) -> FitsResult<()> {
        self.freq().write(fptr, hdu, 1)?;
        self.time().write(fptr, hdu, 2)
    }
}

/// Autocorrelation metrics.
///
/// For each tile with an autocorrelation in the selection, and each
/// polarisation:
/// - the bandpass: the amplitude of the autocorrelation, averaged over the
///   unflagged timesteps;
/// - the autocorrelation amplitude at each timestep, with the bandpass
///   subtracted;
/// - the coefficients of a second order polynomial fit to the bandpass;
/// - the delay transform of the bandpass between 100 and 3000 ns, where
///   reflections in the cables show up as peaks. See
///   [`mod@crate::delay_transform`].
///
/// Flagged channels and timesteps are `NaN`. The flags of individual
/// visibilities are not used.
///
/// The polarisation axes are in the order of the FITS Stokes axis:
/// `XX, YY, XY, YX`. The antenna axes are in order of increasing antenna
/// index.
pub struct AutoMetrics {
    /// The autocorrelation amplitude with its mean over time (the bandpass)
    /// subtracted `(antenna, polarisation, time, frequency)`
    pub auto_sub_aptf: Array4<f32>,
    /// The bandpass: the mean over time of the autocorrelation amplitude
    /// `(antenna, frequency, polarisation)`
    pub auto_spectrum_afp: Array3<f32>,
    /// Coefficients of the polynomial fit to the bandpass, lowest order first,
    /// as a function of frequency in Hz `(antenna, polarisation, order)`. `NaN`
    /// where the fit failed.
    pub auto_coeffs_apo: Array3<f32>,
    /// Delay transform of the bandpass `(antenna, delay, polarisation)`
    pub auto_delay_afp: Array3<f32>,
    /// Tile names of the selected antennas
    pub antenna_names: Vec<String>,
    /// Antenna indices of the selected antennas, which is the order of the
    /// antennas in the metafits. Written to FITS as `ANT_ID`.
    pub antenna_ids: Vec<u32>,
    /// Tile IDs of the selected antennas. Written to FITS as `ANT_NUM`.
    pub antenna_nums: Vec<u32>,
    /// Receiver number of each antenna
    pub rx_numbers: Vec<u32>,
    /// Receiver slot number of each antenna
    pub rx_slots: Vec<u32>,
    /// Receiver type of each antenna
    pub rx_types: Vec<String>,
    /// Cable flavour of each antenna
    pub cable_flavours: Vec<String>,
    /// Whether each antenna has a whitening filter
    pub whitening_filters: Vec<bool>,
    /// Position of each antenna relative to the centre of the array \[m\].
    /// Written to FITS as `EAST`, `NORTH` and `HEIGHT`.
    pub antenna_positions: Vec<ENH>,

    /// Centre frequency of the first channel \[Hz\]
    pub start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    pub channel_width_hz: f64,
    /// GPS time of the first timestep \[s\]
    pub start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    pub integration_time_s: f64,
    /// Delay of the first bin of `auto_delay_afp` \[ns\]
    pub start_delay_ns: f64,
    /// Width of the bins of `auto_delay_afp` \[ns\]
    pub delay_res_ns: f64,
}

/// Computes [`AutoMetrics`] from visibilities which arrive in chunks of
/// timesteps, so that the visibilities of the whole observation don't need to
/// be in memory at once.
///
/// The result is the same as [`AutoMetrics::new_from_metadata`] on all of the
/// visibilities, however they are chunked.
///
/// The amplitude of every autocorrelation is kept, which is `4 * tiles *
/// timesteps * channels` 32-bit floats.
pub struct AutoMetricsBuilder {
    metadata: MetricsContext,
    timestep_flags: Vec<bool>,
    chan_flags: Vec<bool>,
    /// The position in `metadata.antennas`, and the baseline index of each
    /// autocorrelation, in order of increasing antenna index.
    autos: Vec<(usize, usize)>,
    /// The autocorrelation amplitude `(antenna, polarisation, time,
    /// frequency)`, `NaN` where it is flagged or has not been pushed.
    auto_amp_aptf: Array4<f32>,
    num_timesteps_pushed: usize,
}

impl AutoMetricsBuilder {
    /// Start computing [`AutoMetrics`] for all of the timesteps, channels and
    /// baselines in `metadata`.
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - all timesteps or all channels are flagged,
    /// - the lengths of the flags don't match `metadata`.
    pub fn new(metadata: &MetricsContext, timestep_flags: &[bool], chan_flags: &[bool]) -> Self {
        assert!(
            timestep_flags.iter().any(|&t| !t),
            "all timesteps are flagged"
        );
        assert!(chan_flags.iter().any(|&c| !c), "all channels are flagged");

        let sel_ants = selected_antennas(&metadata.antenna_pairs);
        let mut autos: Vec<(usize, usize)> = metadata
            .antenna_pairs
            .iter()
            .enumerate()
            .filter(|&(_, &(a, b))| a == b)
            .map(|(bl_idx, &(a, _))| {
                let ant_pos = sel_ants.binary_search(&a).expect("antenna is in a pair");
                (ant_pos, bl_idx)
            })
            .collect();
        autos.sort_unstable();

        Self {
            auto_amp_aptf: Array4::from_elem(
                (autos.len(), 4, timestep_flags.len(), chan_flags.len()),
                f32::NAN,
            ),
            metadata: metadata.clone(),
            timestep_flags: timestep_flags.to_vec(),
            chan_flags: chan_flags.to_vec(),
            autos,
            num_timesteps_pushed: 0,
        }
    }

    /// Add the visibilities `[timestep][channel][baseline]` of the next
    /// timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the channels or baselines don't match the metadata, or if
    /// there are more timesteps than are left.
    pub fn push(&mut self, jones_array_tfb: ArrayView3<Jones<f32>>) {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let first_timestep = self.num_timesteps_pushed;
        assert!(
            first_timestep + num_timesteps <= self.timestep_flags.len(),
            "too many timesteps"
        );
        assert_eq!(num_freqs, self.chan_flags.len());
        assert_eq!(num_baselines, self.metadata.antenna_pairs.len());

        for (a, &(_, bl_idx)) in self.autos.iter().enumerate() {
            for p in 0..4 {
                let fits_p = hyperdrive_to_fits_stokes(p);
                for t in (0..num_timesteps).filter(|&t| !self.timestep_flags[first_timestep + t]) {
                    for f in (0..num_freqs).filter(|&f| !self.chan_flags[f]) {
                        self.auto_amp_aptf[[a, fits_p, first_timestep + t, f]] =
                            jones_array_tfb[[t, f, bl_idx]][p].norm();
                    }
                }
            }
        }
        self.num_timesteps_pushed += num_timesteps;
    }

    /// Compute the metrics.
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - not all of the timesteps have been pushed,
    /// - there are fewer than 2 channels.
    pub fn finish(self) -> AutoMetrics {
        let Self {
            metadata,
            timestep_flags,
            chan_flags,
            autos,
            auto_amp_aptf: mut auto_sub_aptf,
            num_timesteps_pushed,
        } = self;
        let (num_sel_ants, num_pols, num_timesteps, num_freqs) = auto_sub_aptf.dim();
        assert_eq!(
            num_timesteps_pushed, num_timesteps,
            "not all timesteps were pushed"
        );
        let num_unflagged_times = timestep_flags.iter().filter(|&t| !t).count() as f32;

        // the bandpass: the time-mean of the amplitude
        let mut auto_spectrum_afp =
            Array3::<f32>::from_elem((num_sel_ants, num_freqs, num_pols), f32::NAN);
        for ((a, f, p), mean) in auto_spectrum_afp.indexed_iter_mut() {
            if chan_flags[f] {
                continue;
            }
            *mean = 0.0;
            for t in (0..num_timesteps).filter(|&t| !timestep_flags[t]) {
                *mean += auto_sub_aptf[[a, p, t, f]] / num_unflagged_times;
            }
        }

        // subtract the time-mean
        for ((a, p, _, f), amp) in auto_sub_aptf.indexed_iter_mut() {
            *amp -= auto_spectrum_afp[[a, f, p]];
        }

        // fit a polynomial in frequency to each bandpass
        let mut auto_coeffs_apo =
            Array3::<f32>::from_elem((num_sel_ants, num_pols, AUTO_POLY_ORDER + 1), f32::NAN);
        for a in 0..num_sel_ants {
            for p in 0..num_pols {
                let spectrum = auto_spectrum_afp.slice(s![a, .., p]).to_vec();
                if let Some(coeffs) =
                    fit_polynomial(&metadata.fine_chan_freqs_hz, &spectrum, AUTO_POLY_ORDER)
                {
                    auto_coeffs_apo
                        .slice_mut(s![a, p, ..])
                        .assign(&Array1::from(coeffs));
                }
            }
        }

        // delay transform of each bandpass
        let freqs_hz = Array1::from(metadata.fine_chan_freqs_hz.clone());
        let delay_info = calculate_delay_channels(num_freqs, &freqs_hz, &DELAY_TRANSFORM_CONFIG);
        let mut auto_delay_afp =
            Array3::<f32>::zeros((num_sel_ants, delay_info.n_delay_channels, num_pols));
        for p in 0..num_pols {
            let spectrum_af = auto_spectrum_afp.slice(s![.., .., p]).mapv(f64::from);
            let delay_result = delay_transform(&spectrum_af, &freqs_hz, &DELAY_TRANSFORM_CONFIG)
                .expect("delay_transform failed");
            auto_delay_afp
                .slice_mut(s![.., .., p])
                .assign(&delay_result.delay_spectrum.mapv(|x| x as f32));
        }

        let sel_ant_meta = || {
            autos
                .iter()
                .map(|&(ant_pos, _)| &metadata.antennas[ant_pos])
        };

        AutoMetrics {
            auto_sub_aptf,
            auto_spectrum_afp,
            auto_coeffs_apo,
            auto_delay_afp,
            antenna_names: sel_ant_meta().map(|ant| ant.tile_name.clone()).collect(),
            antenna_ids: sel_ant_meta().map(|ant| ant.ant_id).collect(),
            antenna_nums: sel_ant_meta().map(|ant| ant.tile_id).collect(),
            rx_numbers: sel_ant_meta().map(|ant| ant.rec_number).collect(),
            rx_slots: sel_ant_meta().map(|ant| ant.rec_slot_number).collect(),
            rx_types: sel_ant_meta().map(|ant| ant.rec_type.clone()).collect(),
            cable_flavours: sel_ant_meta()
                .map(|ant| ant.cable_flavour.clone())
                .collect(),
            whitening_filters: sel_ant_meta().map(|ant| ant.has_whitening_filter).collect(),
            antenna_positions: sel_ant_meta()
                .map(|ant| ENH {
                    e: ant.east_m,
                    n: ant.north_m,
                    h: ant.height_m,
                })
                .collect(),
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: metadata.channel_width_hz(),
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: metadata.integration_time_s(),
            start_delay_ns: delay_info.start_delay_ns,
            delay_res_ns: delay_info.actual_delay_res_ns,
        }
    }
}

impl AutoMetrics {
    /// Create [`AutoMetrics`] for visibilities that were read with the
    /// selection `chunk_vis_sel` of the observation `corr_ctx`.
    ///
    /// # Panics
    ///
    /// See [`AutoMetrics::new_from_metadata`].
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let (metadata, timestep_flags, chan_flags) =
            selection_metadata(corr_ctx, chunk_vis_sel, flag_ctx);
        Self::new_from_metadata(jones_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create [`AutoMetrics`] from visibilities `[timestep][channel][baseline]`
    /// and their metadata.
    ///
    /// To compute the metrics from visibilities that don't fit in memory, see
    /// [`AutoMetricsBuilder`].
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - all timesteps or all channels are flagged,
    /// - there are fewer than 2 channels,
    /// - the dimensions of the visibilities or flags don't match `metadata`.
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, _) = jones_array_tfb.dim();
        let mut builder = AutoMetricsBuilder::new(
            metadata,
            &timestep_flags[..num_timesteps],
            &chan_flags[..num_freqs],
        );
        builder.push(jones_array_tfb);
        builder.finish()
    }

    /// Write these metrics to HDUs in an open FITS file:
    ///
    /// - `AUTO_SUB_ANT=<tile name>` for each antenna: an image of
    ///   [`AutoMetrics::auto_sub_aptf`] with FITS axes (frequency, time,
    ///   Stokes), and the metadata of the tile in the header.
    /// - `AUTO_POL=<pol>` for each polarisation: an image of
    ///   [`AutoMetrics::auto_spectrum_afp`] with FITS axes (frequency,
    ///   antenna).
    /// - `AUTO_COEFFS`: a table with a row for each antenna, with the metadata
    ///   of the tile and [`AutoMetrics::auto_coeffs_apo`].
    /// - `AUTO_DELAY_POL=<pol>` for each polarisation: an image of
    ///   [`AutoMetrics::auto_delay_afp`] with FITS axes (delay, antenna).
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the FITS file fails.
    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> FitsResult<()> {
        let (num_ants, num_pols, num_times, num_freqs) = self.auto_sub_aptf.dim();
        let num_delays = self.auto_delay_afp.dim().1;
        let whitening_filters: Vec<i32> = self
            .whitening_filters
            .iter()
            .map(|&w| i32::from(w))
            .collect();

        // every per-antenna field is indexed by `a`
        #[allow(clippy::needless_range_loop)]
        for a in 0..num_ants {
            let hdu = write_image_hdu(
                fptr,
                &format!("AUTO_SUB_ANT={}", self.antenna_names[a]),
                ImageType::Float,
                &[num_pols, num_times, num_freqs],
                self.auto_sub_aptf.slice(s![a, .., .., ..]).iter().copied(),
            )?;
            // CDELT is left off the frequency axis, because carta does not
            // display the image properly with it.
            FitsAxis {
                ctype: "FREQ",
                crval: self.start_freq_hz,
                cdelt: None,
                cunit: Some("Hz"),
            }
            .write(fptr, &hdu, 1)?;
            hdu.write_key(fptr, "CHAN_WIDTH", self.channel_width_hz)?;
            FitsAxis {
                ctype: "TIME",
                crval: self.start_time_gps_s,
                cdelt: Some(self.integration_time_s),
                cunit: Some("s"),
            }
            .write(fptr, &hdu, 2)?;
            FitsAxis {
                ctype: "STOKES",
                crval: -5.0, // XX
                cdelt: Some(-1.0),
                cunit: None,
            }
            .write(fptr, &hdu, 3)?;

            hdu.write_key(fptr, "ANTNAME", self.antenna_names[a].clone())?;
            hdu.write_key(fptr, "ANT_ID", self.antenna_ids[a])?;
            hdu.write_key(fptr, "ANT_NUM", self.antenna_nums[a])?;
            hdu.write_key(fptr, "CABLE_FLAVOUR", self.cable_flavours[a].clone())?;
            hdu.write_key(fptr, "WHITENING_FILTER", whitening_filters[a])?;
            hdu.write_key(fptr, "RX_NUMBER", self.rx_numbers[a])?;
            hdu.write_key(fptr, "RX_SLOT", self.rx_slots[a])?;
            hdu.write_key(fptr, "RX_TYPE", self.rx_types[a].clone())?;

            let position = self.antenna_positions[a];
            hdu.write_key(fptr, "NORTH", position.n)?;
            hdu.write_key(fptr, "EAST", position.e)?;
            hdu.write_key(fptr, "HEIGHT", position.h)?;
        }

        for (p, pol_name) in FITS_POL_NAMES.into_iter().enumerate() {
            let hdu = write_image_hdu(
                fptr,
                &format!("AUTO_POL={pol_name}"),
                ImageType::Float,
                &[num_ants, num_freqs],
                self.auto_spectrum_afp.slice(s![.., .., p]).iter().copied(),
            )?;
            // CDELT is left off the frequency axis, because carta does not
            // display the image properly with it.
            FitsAxis {
                ctype: "FREQ",
                crval: self.start_freq_hz,
                cdelt: None,
                cunit: Some("Hz"),
            }
            .write(fptr, &hdu, 1)?;
            hdu.write_key(fptr, "CHAN_WIDTH", self.channel_width_hz)?;
            hdu.write_key(fptr, "CTYPE2", "BASELINE")?;
            hdu.write_key(fptr, "CRVAL2", 0.0f64)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
        }

        let num_coeffs = self.auto_coeffs_apo.dim().2;
        let string_col = |name| {
            ColumnDescription::new(name)
                .with_type(ColumnDataType::String)
                .that_repeats(32)
                .create()
        };
        let scalar_col =
            |name, data_type| ColumnDescription::new(name).with_type(data_type).create();
        let coeffs_col = |name: &str| {
            ColumnDescription::new(name)
                .with_type(ColumnDataType::Double)
                .that_repeats(num_coeffs)
                .create()
        };
        let mut columns = vec![
            string_col("ANT_NAME")?,
            scalar_col("ANT_ID", ColumnDataType::Int)?,
            scalar_col("ANT_NUM", ColumnDataType::Int)?,
            string_col("CABLE_FLAVOUR")?,
            scalar_col("WHITENING_FILTER", ColumnDataType::Int)?,
            scalar_col("RX_NUMBER", ColumnDataType::Int)?,
            scalar_col("RX_SLOT", ColumnDataType::Int)?,
            string_col("RX_TYPE")?,
            scalar_col("NORTH", ColumnDataType::Double)?,
            scalar_col("EAST", ColumnDataType::Double)?,
            scalar_col("HEIGHT", ColumnDataType::Double)?,
        ];
        for pol_name in FITS_POL_NAMES {
            columns.push(coeffs_col(&format!("COEFFS_{pol_name}"))?);
        }

        let hdu = fptr.create_table("AUTO_COEFFS", &columns)?;
        hdu.write_col(fptr, "ANT_NAME", &self.antenna_names)?;
        hdu.write_col(fptr, "ANT_ID", &self.antenna_ids)?;
        hdu.write_col(fptr, "ANT_NUM", &self.antenna_nums)?;
        hdu.write_col(fptr, "CABLE_FLAVOUR", &self.cable_flavours)?;
        hdu.write_col(fptr, "WHITENING_FILTER", &whitening_filters)?;
        hdu.write_col(fptr, "RX_NUMBER", &self.rx_numbers)?;
        hdu.write_col(fptr, "RX_SLOT", &self.rx_slots)?;
        hdu.write_col(fptr, "RX_TYPE", &self.rx_types)?;
        let positions = &self.antenna_positions;
        hdu.write_col(
            fptr,
            "NORTH",
            &positions.iter().map(|p| p.n).collect::<Vec<_>>(),
        )?;
        hdu.write_col(
            fptr,
            "EAST",
            &positions.iter().map(|p| p.e).collect::<Vec<_>>(),
        )?;
        hdu.write_col(
            fptr,
            "HEIGHT",
            &positions.iter().map(|p| p.h).collect::<Vec<_>>(),
        )?;
        for (p, pol_name) in FITS_POL_NAMES.into_iter().enumerate() {
            // one row of coefficients per antenna
            let coeffs: Vec<f64> = self
                .auto_coeffs_apo
                .slice(s![.., p, ..])
                .iter()
                .map(|&c| f64::from(c))
                .collect();
            hdu.write_col(fptr, format!("COEFFS_{pol_name}"), &coeffs)?;
        }
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "POLY_ORD", (num_coeffs - 1) as u32)?;

        for (p, pol_name) in FITS_POL_NAMES.into_iter().enumerate() {
            let hdu = write_image_hdu(
                fptr,
                &format!("AUTO_DELAY_POL={pol_name}"),
                ImageType::Float,
                &[num_ants, num_delays],
                self.auto_delay_afp.slice(s![.., .., p]).iter().copied(),
            )?;
            FitsAxis {
                ctype: "DELAY",
                crval: self.start_delay_ns,
                cdelt: None,
                cunit: Some("ns"),
            }
            .write(fptr, &hdu, 1)?;
            hdu.write_key(fptr, "DELAY_RES", self.delay_res_ns)?;
            hdu.write_key(fptr, "CTYPE2", "BASELINE")?;
            hdu.write_key(fptr, "CRVAL2", 0.0f64)?;
            hdu.write_key(fptr, "CRPIX2", 1.0f64)?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
        }

        Ok(())
    }
}

/// SSINS (Sky-Subtracted Incoherent Noise Spectra)
///
/// The amplitude of the difference between the visibilities at adjacent
/// timesteps, averaged incoherently over the cross-correlation baselines, as a
/// z-score. The sky changes slowly, so it cancels in the difference, leaving
/// noise and anything that changes quickly, like radio-frequency interference.
///
/// It is computed like `SSINS.INS(SS.read(..., diff=True), spectrum_type="cross", order=0)`
/// in [SSINS](https://ssins.readthedocs.io/en/latest/):
/// - only cross-correlations are averaged;
/// - a baseline contributes to (t, f) only if neither of the visibilities being
///   differenced is flagged, so `N_bl` is counted per (t, f);
/// - the time mean is weighted by `N_bl`.
///
/// The product is a 3D z-score array with dimensions
/// `(num_timesteps - 1, num_freqs, num_pols=4)`, where
///
/// ```text
/// C = 4/pi - 1
/// zscore = sqrt(N_bl(t, f) / C) * (mean_amp_tfp / mean_amp_fp - 1)
/// ```
///
/// Cells with no unflagged baselines are `NaN`. The polarisations are in the
/// order of a [`Jones`] matrix: `XX, XY, YX, YY`.
#[allow(clippy::upper_case_acronyms)]
pub struct SSINS {
    /// z-score of the baseline-averaged amplitude of the time difference
    /// `(times-1, frequencies, polarisations)`
    pub zscore: Array3<f32>,
    /// The mean over time of the baseline-averaged amplitude of the time
    /// difference, weighted by the number of baselines
    /// `(frequencies, polarisations)`
    pub diff_mean_amp_fp: Array2<f32>,
    /// Flags for the z-scores `(times-1, frequencies)`. These are all `false`,
    /// unless they have been set by `SSINS::flag` (`aoflagger` feature).
    pub flag_array: Array2<bool>,
    /// number of unflagged cross-correlation baselines averaged in each (time-1, frequency)
    pub num_baselines_tf: Array2<u32>,
    /// number of cross-correlation baselines selected
    pub num_baselines: usize,
    /// GPS time halfway between the first two timesteps \[s\]
    pub start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    pub integration_time_s: f64,
    /// Centre frequency of the first channel \[Hz\]
    pub start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    pub channel_width_hz: f64,
}

/// Computes [`SSINS`] from visibilities which arrive in chunks of timesteps,
/// so that the visibilities of the whole observation don't need to be in
/// memory at once.
///
/// The result is the same as [`SSINS::new_from_metadata`] on all of the
/// visibilities, however they are chunked: the last timestep of each chunk is
/// kept, to difference with the first timestep of the next.
pub struct SSINSBuilder {
    metadata: MetricsContext,
    timestep_flags: Vec<bool>,
    chan_flags: Vec<bool>,
    /// Baseline indices of the cross-correlations. SSINS's
    /// `spectrum_type="cross"` leaves out the autos.
    cross_bl_idxs: Vec<usize>,
    /// Sum over unflagged cross baselines of `|V(t+1) - V(t)|`
    /// `(times-1, frequencies, polarisations)`
    diff_sum_amp_tfp: Array3<f64>,
    /// How many baselines are in each sum `(times-1, frequencies)`
    num_baselines_tf: Array2<u32>,
    /// The visibilities and flags `[channel][baseline]` of the last timestep
    /// that was pushed, if there are more to come.
    last_timestep: Option<(Array2<Jones<f32>>, Array2<bool>)>,
    num_timesteps_pushed: usize,
}

impl SSINSBuilder {
    /// Start computing [`SSINS`] for all of the timesteps, channels and
    /// baselines in `metadata`.
    ///
    /// # Panics
    ///
    /// Panics if there are fewer than 2 timesteps.
    pub fn new(metadata: &MetricsContext, timestep_flags: &[bool], chan_flags: &[bool]) -> Self {
        let num_timesteps = timestep_flags.len();
        let num_freqs = chan_flags.len();
        assert!(num_timesteps >= 2, "SSINS requires at least 2 timesteps");

        let cross_bl_idxs: Vec<usize> = metadata
            .antenna_pairs
            .iter()
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(i, _)| i)
            .collect();

        Self {
            metadata: metadata.clone(),
            timestep_flags: timestep_flags.to_vec(),
            chan_flags: chan_flags.to_vec(),
            cross_bl_idxs,
            diff_sum_amp_tfp: Array3::zeros((num_timesteps - 1, num_freqs, 4)),
            num_baselines_tf: Array2::zeros((num_timesteps - 1, num_freqs)),
            last_timestep: None,
            num_timesteps_pushed: 0,
        }
    }

    /// Add the visibilities `[timestep][channel][baseline]` and flags of the
    /// next timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the visibilities and flags don't have the same dimensions,
    /// the channels or baselines don't match the metadata, or if there are
    /// more timesteps than are left.
    pub fn push(
        &mut self,
        jones_array_tfb: ArrayView3<Jones<f32>>,
        flag_array_tfb: ArrayView3<bool>,
    ) {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        assert_eq!(jones_array_tfb.dim(), flag_array_tfb.dim());
        let first_timestep = self.num_timesteps_pushed;
        let num_timesteps_total = self.timestep_flags.len();
        assert!(
            first_timestep + num_timesteps <= num_timesteps_total,
            "too many timesteps"
        );
        assert_eq!(num_freqs, self.chan_flags.len());
        assert_eq!(num_baselines, self.metadata.antenna_pairs.len());
        if num_timesteps == 0 {
            return;
        }

        let Self {
            timestep_flags,
            chan_flags,
            cross_bl_idxs,
            diff_sum_amp_tfp,
            num_baselines_tf,
            ..
        } = self;

        // The difference between timestep `t` and `t + 1` of the observation
        let diff_timesteps =
            |t: usize,
             (jones_fb, flags_fb): (ArrayView2<Jones<f32>>, ArrayView2<bool>),
             (next_jones_fb, next_flags_fb): (ArrayView2<Jones<f32>>, ArrayView2<bool>),
             mut sum_fp: ArrayViewMut2<f64>,
             mut count_f: ArrayViewMut1<u32>| {
                if timestep_flags[t] || timestep_flags[t + 1] {
                    return;
                }
                for f in 0..num_freqs {
                    if chan_flags[f] {
                        continue;
                    }
                    for &b in cross_bl_idxs.iter() {
                        if flags_fb[[f, b]] || next_flags_fb[[f, b]] {
                            continue;
                        }
                        let diff = next_jones_fb[[f, b]] - jones_fb[[f, b]];
                        for p in 0..4 {
                            sum_fp[[f, p]] += f64::from(diff[p].norm());
                        }
                        count_f[f] += 1;
                    }
                }
            };
        let timestep = |t: usize| {
            (
                jones_array_tfb.index_axis(Axis(0), t),
                flag_array_tfb.index_axis(Axis(0), t),
            )
        };

        // the last timestep of the previous chunk, with the first of this one
        if let Some((last_jones_fb, last_flags_fb)) = &self.last_timestep {
            diff_timesteps(
                first_timestep - 1,
                (last_jones_fb.view(), last_flags_fb.view()),
                timestep(0),
                diff_sum_amp_tfp.index_axis_mut(Axis(0), first_timestep - 1),
                num_baselines_tf.index_axis_mut(Axis(0), first_timestep - 1),
            );
        }

        // the timesteps within this chunk
        let diff_range = first_timestep..first_timestep + num_timesteps - 1;
        diff_sum_amp_tfp
            .slice_mut(s![diff_range.clone(), .., ..])
            .outer_iter_mut()
            .into_par_iter()
            .zip(
                num_baselines_tf
                    .slice_mut(s![diff_range, ..])
                    .outer_iter_mut(),
            )
            .enumerate()
            .for_each(|(t, (sum_fp, count_f))| {
                diff_timesteps(
                    first_timestep + t,
                    timestep(t),
                    timestep(t + 1),
                    sum_fp,
                    count_f,
                );
            });

        self.num_timesteps_pushed += num_timesteps;
        self.last_timestep = (self.num_timesteps_pushed < num_timesteps_total).then(|| {
            let (jones_fb, flags_fb) = timestep(num_timesteps - 1);
            (jones_fb.to_owned(), flags_fb.to_owned())
        });
    }

    /// Compute the metrics.
    ///
    /// # Panics
    ///
    /// Panics if not all of the timesteps have been pushed.
    pub fn finish(self) -> SSINS {
        let Self {
            metadata,
            cross_bl_idxs,
            diff_sum_amp_tfp,
            num_baselines_tf,
            num_timesteps_pushed,
            timestep_flags,
            ..
        } = self;
        let (num_diffs, num_freqs, _) = diff_sum_amp_tfp.dim();
        assert_eq!(
            num_timesteps_pushed,
            timestep_flags.len(),
            "not all timesteps were pushed"
        );

        // time mean weighted by the number of baselines, like SSINS's np.ma.average.
        let mut diff_mean_amp_fp = Array2::<f32>::from_elem((num_freqs, 4), f32::NAN);
        for f in 0..num_freqs {
            let count: f64 = num_baselines_tf
                .column(f)
                .iter()
                .map(|&n| f64::from(n))
                .sum();
            if count > 0.0 {
                for p in 0..4 {
                    diff_mean_amp_fp[[f, p]] =
                        (diff_sum_amp_tfp.slice(s![.., f, p]).sum() / count) as f32;
                }
            }
        }

        let c = 4.0 / std::f64::consts::PI - 1.0;
        let mut zscore = Array3::<f32>::from_elem((num_diffs, num_freqs, 4), f32::NAN);
        for t in 0..num_diffs {
            for f in 0..num_freqs {
                let n = num_baselines_tf[[t, f]];
                if n == 0 {
                    continue;
                }
                let n = f64::from(n);
                for p in 0..4 {
                    let mean_amp = diff_sum_amp_tfp[[t, f, p]] / n;
                    zscore[[t, f, p]] = ((n / c).sqrt()
                        * (mean_amp / f64::from(diff_mean_amp_fp[[f, p]]) - 1.0))
                        as f32;
                }
            }
        }

        SSINS {
            zscore,
            diff_mean_amp_fp,
            flag_array: Array2::<bool>::default((num_diffs, num_freqs)),
            num_baselines_tf,
            num_baselines: cross_bl_idxs.len(),
            // each difference is halfway between two timesteps
            start_time_gps_s: (metadata.timestamps_s[0] + metadata.timestamps_s[1]) / 2.0,
            integration_time_s: metadata.integration_time_s(),
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: metadata.channel_width_hz(),
        }
    }
}

impl SSINS {
    /// Create [`SSINS`] for visibilities that were read with the selection
    /// `chunk_vis_sel` of the observation `corr_ctx`.
    ///
    /// # Panics
    ///
    /// See [`SSINS::new_from_metadata`].
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        flag_array_tfb: ArrayView3<bool>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let (metadata, timestep_flags, chan_flags) =
            selection_metadata(corr_ctx, chunk_vis_sel, flag_ctx);
        Self::new_from_metadata(
            jones_array_tfb,
            flag_array_tfb,
            &metadata,
            &timestep_flags,
            &chan_flags,
        )
    }

    /// Create [`SSINS`] from visibilities `[timestep][channel][baseline]`,
    /// their flags, and their metadata.
    ///
    /// To compute the metrics from visibilities that don't fit in memory, see
    /// [`SSINSBuilder`].
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - there are fewer than 2 timesteps,
    /// - the visibilities and flags don't have the same dimensions,
    /// - the dimensions of the visibilities or flags don't match `metadata`.
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        flag_array_tfb: ArrayView3<bool>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, _) = jones_array_tfb.dim();
        let mut builder = SSINSBuilder::new(
            metadata,
            &timestep_flags[..num_timesteps],
            &chan_flags[..num_freqs],
        );
        builder.push(jones_array_tfb, flag_array_tfb);
        builder.finish()
    }

    /// Flag the z-scores with an `AOFlagger` strategy (the generic "minimal"
    /// strategy if none is given), writing the result into `flag_array`.
    #[cfg(feature = "aoflagger")]
    pub fn flag(&mut self, strategy_filename: Option<String>) {
        use aoflagger_sys::cxx_aoflagger_new;

        use crate::flags::{amps_tfp_to_imageset, flag_baseline_view_to_flagmask};

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let imgset = amps_tfp_to_imageset(&aoflagger, self.zscore.view());
        let strategy_filename = strategy_filename
            .unwrap_or_else(|| aoflagger.FindStrategyFileGeneric(&String::from("minimal")));
        let flag_strategy = aoflagger.LoadStrategyFile(&strategy_filename);

        let flagmask = flag_baseline_view_to_flagmask(&aoflagger, self.flag_array.view());
        let new_flagmask = flag_strategy.RunExisting(&imgset, &flagmask);
        let flag_buf = new_flagmask.Buffer();
        let stride = new_flagmask.HorizontalStride();
        for ((img_timestep_idx, img_chan_idx), flag) in self.flag_array.indexed_iter_mut() {
            *flag = flag_buf[img_chan_idx * stride + img_timestep_idx];
        }
    }

    /// Write these metrics to HDUs in an open FITS file:
    ///
    /// - `SSINS_POL=<pol>` for each polarisation: an image of [`SSINS::zscore`]
    ///   with FITS axes (frequency, time).
    /// - `SSINS_FLAGS`: an image of [`SSINS::flag_array`] (1 is flagged) with
    ///   FITS axes (frequency, time).
    /// - `SSINS_N_BL`: an image of [`SSINS::num_baselines_tf`] with FITS axes
    ///   (frequency, time).
    /// - `SSINS_DIFF_MEAN_AMP_FP`: an image of [`SSINS::diff_mean_amp_fp`] with
    ///   FITS axes (polarisation, frequency).
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the FITS file fails.
    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> FitsResult<()> {
        let (num_times, num_freqs, _) = self.zscore.dim();
        let axes = TimeFreqAxes {
            start_time_gps_s: self.start_time_gps_s,
            integration_time_s: self.integration_time_s,
            start_freq_hz: self.start_freq_hz,
            channel_width_hz: self.channel_width_hz,
        };

        for (p, pol_name) in JONES_POL_NAMES.into_iter().enumerate() {
            let hdu = write_image_hdu(
                fptr,
                &format!("SSINS_POL={pol_name}"),
                ImageType::Float,
                &[num_times, num_freqs],
                self.zscore.slice(s![.., .., p]).iter().copied(),
            )?;
            axes.write(fptr, &hdu)?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_BL", self.num_baselines as u32)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "SSINS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;
            hdu.write_key(
                fptr,
                "COMMENT",
                "SSINS (Sky-Subtracted Incoherent Noise Spectra)",
            )?;
            hdu.write_key(
                fptr,
                "COMMENT",
                "Incoherently averaged difference between visibilities in time",
            )?;
            hdu.write_key(
                fptr,
                "COMMENT",
                "One image per polarization (XX, XY, YX, YY)",
            )?;
        }

        let hdu = write_image_hdu(
            fptr,
            "SSINS_FLAGS",
            ImageType::UnsignedByte,
            &[num_times, num_freqs],
            self.flag_array.iter().map(|&b| u8::from(b)),
        )?;
        axes.write(fptr, &hdu)?;
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "SSINS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        // the number of unflagged cross baselines in each cell, which SSINS
        // keeps as its weights_array.
        let hdu = write_image_hdu(
            fptr,
            "SSINS_N_BL",
            ImageType::Long,
            &[num_times, num_freqs],
            self.num_baselines_tf.iter().map(|&n| n as i32),
        )?;
        axes.write(fptr, &hdu)?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;

        let hdu = write_image_hdu(
            fptr,
            "SSINS_DIFF_MEAN_AMP_FP",
            ImageType::Float,
            &[num_freqs, 4],
            self.diff_mean_amp_fp.iter().copied(),
        )?;
        // the polarisations are FITS axis 1, and frequency is FITS axis 2.
        axes.freq().write(fptr, &hdu, 2)?;

        Ok(())
    }
}

/// EAVILS (Expected Amplitude of `VisibILities` Spectra)
///
/// The amplitude of the visibilities, averaged incoherently over all selected
/// baselines (including autocorrelations), as a z-score relative to its mean
/// over time. This is like [`SSINS`], but without differencing the
/// visibilities in time, so it is also sensitive to things that don't change
/// between timesteps.
///
/// With `amp = |V|` for visibilities with dimensions (time, baseline,
/// frequency, polarisation), it is equivalent to
///
/// ```python
/// blmean = np.mean(amp, axis=1)
/// stdv = np.sqrt(np.mean(np.var(amp, axis=0, ddof=1), axis=0))
/// zscore = (blmean - np.mean(blmean, axis=0)) / stdv * np.sqrt(Nbls)
/// ```
///
/// where the means and variances over time only include unflagged timesteps,
/// and the variance is the sample variance (it divides by `N - 1`). With
/// fewer than two unflagged timesteps there is no variance, and the z-scores
/// are not finite.
/// Flagged channels and timesteps are `NaN`. The flags of individual
/// visibilities are not used.
///
/// The polarisations are in the order of a [`Jones`] matrix: `XX, XY, YX, YY`.
#[allow(clippy::upper_case_acronyms)]
pub struct EAVILS {
    /// `sqrt(Nbls) * (mean_amp_tfp - mean_amp_fp) / sqrt_mean_var_amp_fp`
    /// `(times, frequencies, polarisations)`
    pub zscore: Array3<f32>,
    /// The mean over time of the baseline-averaged amplitude
    /// `(frequencies, polarisations)`
    pub mean_amp_fp: Array2<f32>,
    /// The square root of the baseline-averaged sample variance over time of
    /// the amplitude `(frequencies, polarisations)`
    pub sqrt_mean_var_amp_fp: Array2<f32>,
    /// Flags for the z-scores `(times, frequencies)`. These are all `false`.
    pub flag_array: Array2<bool>,
    /// GPS time of the first timestep \[s\]
    pub start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    pub integration_time_s: f64,
    /// Centre frequency of the first channel \[Hz\]
    pub start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    pub channel_width_hz: f64,
}

/// Computes [`EAVILS`] from visibilities which arrive in chunks of timesteps,
/// so that the visibilities of the whole observation don't need to be in
/// memory at once.
///
/// The result does not depend on how the visibilities are chunked: the
/// variance over time of each baseline is accumulated one timestep at a time
/// with Welford's algorithm.
///
/// The mean and variance of every polarisation, baseline and channel are kept
/// as 64-bit floats, because accumulating them in 32 bits loses precision.
pub struct EAVILSBuilder {
    metadata: MetricsContext,
    timestep_flags: Vec<bool>,
    chan_flags: Vec<bool>,
    /// Sum over baselines of the amplitude `(times, frequencies,
    /// polarisations)`, `NaN` where it is flagged.
    sum_amp_tfp: Array3<f32>,
    /// Mean over the unflagged timesteps so far of the amplitude
    /// `(frequencies, baselines, polarisations)`
    mean_amp_fbp: Array3<f64>,
    /// Sum over the unflagged timesteps so far of the squared difference of
    /// the amplitude from its mean `(frequencies, baselines, polarisations)`
    sum_sq_dev_amp_fbp: Array3<f64>,
    num_unflagged_timesteps_pushed: usize,
    num_timesteps_pushed: usize,
}

impl EAVILSBuilder {
    /// Start computing [`EAVILS`] for all of the timesteps, channels and
    /// baselines in `metadata`.
    pub fn new(metadata: &MetricsContext, timestep_flags: &[bool], chan_flags: &[bool]) -> Self {
        let num_timesteps = timestep_flags.len();
        let num_freqs = chan_flags.len();
        let num_baselines = metadata.antenna_pairs.len();
        Self {
            metadata: metadata.clone(),
            timestep_flags: timestep_flags.to_vec(),
            chan_flags: chan_flags.to_vec(),
            sum_amp_tfp: Array3::zeros((num_timesteps, num_freqs, 4)),
            mean_amp_fbp: Array3::zeros((num_freqs, num_baselines, 4)),
            sum_sq_dev_amp_fbp: Array3::zeros((num_freqs, num_baselines, 4)),
            num_unflagged_timesteps_pushed: 0,
            num_timesteps_pushed: 0,
        }
    }

    /// Add the visibilities `[timestep][channel][baseline]` of the next
    /// timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the channels or baselines don't match the metadata, or if
    /// there are more timesteps than are left.
    pub fn push(&mut self, jones_array_tfb: ArrayView3<Jones<f32>>) {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let first_timestep = self.num_timesteps_pushed;
        assert!(
            first_timestep + num_timesteps <= self.timestep_flags.len(),
            "too many timesteps"
        );
        assert_eq!(num_freqs, self.chan_flags.len());
        assert_eq!(num_baselines, self.metadata.antenna_pairs.len());

        for t in 0..num_timesteps {
            let mut sum_amp_fp = self.sum_amp_tfp.index_axis_mut(Axis(0), first_timestep + t);
            if self.timestep_flags[first_timestep + t] {
                sum_amp_fp.fill(f32::NAN);
                continue;
            }
            self.num_unflagged_timesteps_pushed += 1;
            let n = self.num_unflagged_timesteps_pushed as f64;
            for f in 0..num_freqs {
                if self.chan_flags[f] {
                    sum_amp_fp.row_mut(f).fill(f32::NAN);
                    continue;
                }
                for b in 0..num_baselines {
                    for p in 0..4 {
                        let amp = jones_array_tfb[[t, f, b]][p].norm();
                        sum_amp_fp[[f, p]] += amp;
                        // Welford's algorithm
                        let mean = &mut self.mean_amp_fbp[[f, b, p]];
                        let amp = f64::from(amp);
                        let dev = amp - *mean;
                        *mean += dev / n;
                        self.sum_sq_dev_amp_fbp[[f, b, p]] += dev * (amp - *mean);
                    }
                }
            }
        }
        self.num_timesteps_pushed += num_timesteps;
    }

    /// Compute the metrics.
    ///
    /// # Panics
    ///
    /// Panics if not all of the timesteps have been pushed.
    pub fn finish(self) -> EAVILS {
        let Self {
            metadata,
            timestep_flags,
            chan_flags,
            sum_amp_tfp: mut mean_amp_tfp,
            sum_sq_dev_amp_fbp: mut var_amp_fbp,
            num_unflagged_timesteps_pushed: num_unflagged_timesteps,
            num_timesteps_pushed,
            ..
        } = self;
        let (num_timesteps, num_freqs, _) = mean_amp_tfp.dim();
        let num_baselines = metadata.antenna_pairs.len();
        assert_eq!(
            num_timesteps_pushed, num_timesteps,
            "not all timesteps were pushed"
        );

        mean_amp_tfp /= num_baselines as f32;

        // the mean over time of the baseline-mean
        let mut mean_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        for t in (0..num_timesteps).filter(|&t| !timestep_flags[t]) {
            for f in 0..num_freqs {
                for p in 0..4 {
                    mean_amp_fp[[f, p]] += mean_amp_tfp[[t, f, p]];
                }
            }
        }
        mean_amp_fp /= num_unflagged_timesteps as f32;

        // the time-variance of the amplitude for each freq, bl, pol
        // (the sample variance, like numpy's `ddof=1`)
        if num_unflagged_timesteps > 1 {
            var_amp_fbp /= (num_unflagged_timesteps - 1) as f64;
        }

        // the square root of the baseline-mean of the variance
        let mut sqrt_mean_var_amp_fp = Array2::<f32>::zeros((num_freqs, 4));
        for f in 0..num_freqs {
            for p in 0..4 {
                if chan_flags[f] {
                    mean_amp_fp[[f, p]] = f32::NAN;
                    sqrt_mean_var_amp_fp[[f, p]] = f32::NAN;
                    continue;
                }
                for b in 0..num_baselines {
                    sqrt_mean_var_amp_fp[[f, p]] += var_amp_fbp[[f, b, p]] as f32;
                }
            }
        }
        sqrt_mean_var_amp_fp /= num_baselines as f32;
        sqrt_mean_var_amp_fp = sqrt_mean_var_amp_fp.sqrt();

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

        EAVILS {
            zscore,
            mean_amp_fp,
            sqrt_mean_var_amp_fp,
            flag_array: Array2::<bool>::default((num_timesteps, num_freqs)),
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: metadata.integration_time_s(),
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: metadata.channel_width_hz(),
        }
    }
}

impl EAVILS {
    /// Create [`EAVILS`] for visibilities that were read with the selection
    /// `chunk_vis_sel` of the observation `corr_ctx`.
    ///
    /// # Panics
    ///
    /// See [`EAVILS::new_from_metadata`].
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let (metadata, timestep_flags, chan_flags) =
            selection_metadata(corr_ctx, chunk_vis_sel, flag_ctx);
        Self::new_from_metadata(jones_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create [`EAVILS`] from visibilities `[timestep][channel][baseline]` and
    /// their metadata.
    ///
    /// To compute the metrics from visibilities that don't fit in memory, see
    /// [`EAVILSBuilder`].
    ///
    /// # Panics
    ///
    /// Panics if the dimensions of the visibilities or flags don't match
    /// `metadata`.
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, _) = jones_array_tfb.dim();
        let mut builder = EAVILSBuilder::new(
            metadata,
            &timestep_flags[..num_timesteps],
            &chan_flags[..num_freqs],
        );
        builder.push(jones_array_tfb);
        builder.finish()
    }

    /// Write these metrics to HDUs in an open FITS file:
    ///
    /// - `EAVILS_POL=<pol>` for each polarisation: an image of
    ///   [`EAVILS::zscore`] with FITS axes (frequency, time).
    /// - `EAVILS_FLAGS`: an image of [`EAVILS::flag_array`] (1 is flagged) with
    ///   FITS axes (frequency, time).
    /// - `EAVILS_MEAN_AMP_FP`: an image of [`EAVILS::mean_amp_fp`] with FITS
    ///   axes (polarisation, frequency).
    /// - `EAVILS_SQRT_MEAN_VAR_AMP_FP`: an image of
    ///   [`EAVILS::sqrt_mean_var_amp_fp`] with FITS axes (polarisation,
    ///   frequency).
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the FITS file fails.
    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> FitsResult<()> {
        const DESCRIPTION: &str = "EAVILS (Expected Amplitude of VisibILities Spectra)";
        let (num_times, num_freqs, _) = self.zscore.dim();
        let axes = TimeFreqAxes {
            start_time_gps_s: self.start_time_gps_s,
            integration_time_s: self.integration_time_s,
            start_freq_hz: self.start_freq_hz,
            channel_width_hz: self.channel_width_hz,
        };

        for (p, pol_name) in JONES_POL_NAMES.into_iter().enumerate() {
            let hdu = write_image_hdu(
                fptr,
                &format!("EAVILS_POL={pol_name}"),
                ImageType::Float,
                &[num_times, num_freqs],
                self.zscore.slice(s![.., .., p]).iter().copied(),
            )?;
            axes.write(fptr, &hdu)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "COMMENT", DESCRIPTION)?;
        }

        let hdu = write_image_hdu(
            fptr,
            "EAVILS_FLAGS",
            ImageType::UnsignedByte,
            &[num_times, num_freqs],
            self.flag_array.iter().map(|&b| u8::from(b)),
        )?;
        axes.write(fptr, &hdu)?;
        hdu.write_key(fptr, "TELESCOP", "MWA")?;
        hdu.write_key(fptr, "INSTRUME", "EAVILS")?;
        hdu.write_key(fptr, "ORIGIN", "Birli")?;
        hdu.write_key(fptr, "COMMENT", DESCRIPTION)?;

        // the polarisations are FITS axis 1, and frequency is FITS axis 2.
        for (extname, array_fp) in [
            ("EAVILS_MEAN_AMP_FP", &self.mean_amp_fp),
            ("EAVILS_SQRT_MEAN_VAR_AMP_FP", &self.sqrt_mean_var_amp_fp),
        ] {
            let hdu = write_image_hdu(
                fptr,
                extname,
                ImageType::Float,
                &[num_freqs, 4],
                array_fp.iter().copied(),
            )?;
            axes.freq().write(fptr, &hdu, 2)?;
        }

        Ok(())
    }
}

/// Flag occupancy metrics: the fraction of the selected cross-correlation
/// baselines that are flagged at each time and frequency, for example by
/// `AOFlagger`.
#[cfg(feature = "aoflagger")]
pub struct AOFlagMetrics {
    /// Fraction of cross-correlations which are flagged `(times, frequencies)`.
    /// `NaN` for flagged timesteps and channels, and if there are no
    /// cross-correlations.
    pub occupancy_tf: Array2<f64>,
    /// GPS time of the first timestep \[s\]
    pub start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    pub integration_time_s: f64,
    /// Centre frequency of the first channel \[Hz\]
    pub start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    pub channel_width_hz: f64,
}

/// Computes [`AOFlagMetrics`] from flags which arrive in chunks of timesteps,
/// so that the flags of the whole observation don't need to be in memory at
/// once.
///
/// The result is the same as [`AOFlagMetrics::new_from_metadata`] on all of
/// the flags, however they are chunked.
#[cfg(feature = "aoflagger")]
pub struct AOFlagMetricsBuilder {
    metadata: MetricsContext,
    timestep_flags: Vec<bool>,
    chan_flags: Vec<bool>,
    occupancy_tf: Array2<f64>,
    num_timesteps_pushed: usize,
}

#[cfg(feature = "aoflagger")]
impl AOFlagMetricsBuilder {
    /// Start computing [`AOFlagMetrics`] for all of the timesteps, channels
    /// and baselines in `metadata`.
    pub fn new(metadata: &MetricsContext, timestep_flags: &[bool], chan_flags: &[bool]) -> Self {
        Self {
            metadata: metadata.clone(),
            timestep_flags: timestep_flags.to_vec(),
            chan_flags: chan_flags.to_vec(),
            occupancy_tf: Array2::from_elem((timestep_flags.len(), chan_flags.len()), f64::NAN),
            num_timesteps_pushed: 0,
        }
    }

    /// Add the flags `[timestep][channel][baseline]` of the next timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the channels or baselines don't match the metadata, or if
    /// there are more timesteps than are left.
    pub fn push(&mut self, flag_array_tfb: ArrayView3<bool>) {
        let (num_timesteps, num_freqs, num_baselines) = flag_array_tfb.dim();
        let timestep_range = self.num_timesteps_pushed..self.num_timesteps_pushed + num_timesteps;
        assert!(
            timestep_range.end <= self.timestep_flags.len(),
            "too many timesteps"
        );
        assert_eq!(num_freqs, self.chan_flags.len());
        assert_eq!(num_baselines, self.metadata.antenna_pairs.len());

        self.occupancy_tf
            .slice_mut(s![timestep_range.clone(), ..])
            .assign(&AOFlagMetrics::compute_cross_occupancy(
                flag_array_tfb,
                &self.timestep_flags[timestep_range],
                &self.chan_flags,
                &self.metadata.antenna_pairs,
            ));
        self.num_timesteps_pushed += num_timesteps;
    }

    /// Compute the metrics.
    ///
    /// # Panics
    ///
    /// Panics if not all of the timesteps have been pushed.
    pub fn finish(self) -> AOFlagMetrics {
        assert_eq!(
            self.num_timesteps_pushed,
            self.timestep_flags.len(),
            "not all timesteps were pushed"
        );
        let metadata = &self.metadata;
        AOFlagMetrics {
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: metadata.integration_time_s(),
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: metadata.channel_width_hz(),
            occupancy_tf: self.occupancy_tf,
        }
    }
}

#[cfg(feature = "aoflagger")]
impl AOFlagMetrics {
    /// Compute occupancy across baselines excluding autocorrelations
    ///
    /// `occupancy[t,f]` = fraction of cross-correlation baselines flagged at `(t,f)`
    fn compute_cross_occupancy(
        flag_array_tfb: ArrayView3<bool>,
        timestep_flags: &[bool],
        chan_flags: &[bool],
        ant_pairs: &[(usize, usize)],
    ) -> Array2<f64> {
        let (num_timesteps, num_freqs, _) = flag_array_tfb.dim();

        let cross_baseline_indices: Vec<usize> = ant_pairs
            .iter()
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(i, _)| i)
            .collect();
        let num_cross_baselines = cross_baseline_indices.len();

        Array2::from_shape_fn((num_timesteps, num_freqs), |(t, f)| {
            if timestep_flags[t] || chan_flags[f] || num_cross_baselines == 0 {
                return f64::NAN;
            }
            let num_flagged = cross_baseline_indices
                .iter()
                .filter(|&&b| flag_array_tfb[[t, f, b]])
                .count();
            num_flagged as f64 / num_cross_baselines as f64
        })
    }

    /// Create [`AOFlagMetrics`] for flags of the selection `chunk_vis_sel` of
    /// the observation `corr_ctx`.
    pub fn new(
        flag_array_tfb: ArrayView3<bool>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
    ) -> Self {
        let (metadata, timestep_flags, chan_flags) =
            selection_metadata(corr_ctx, chunk_vis_sel, flag_ctx);
        Self::new_from_metadata(flag_array_tfb, &metadata, &timestep_flags, &chan_flags)
    }

    /// Create [`AOFlagMetrics`] from flags `[timestep][channel][baseline]` and
    /// their metadata.
    ///
    /// To compute the metrics from flags that don't fit in memory, see
    /// [`AOFlagMetricsBuilder`].
    ///
    /// # Panics
    ///
    /// Panics if the dimensions of the flags don't match `metadata`.
    pub fn new_from_metadata(
        flag_array_tfb: ArrayView3<bool>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
    ) -> Self {
        let (num_timesteps, num_freqs, _) = flag_array_tfb.dim();
        let mut builder = AOFlagMetricsBuilder::new(
            metadata,
            &timestep_flags[..num_timesteps],
            &chan_flags[..num_freqs],
        );
        builder.push(flag_array_tfb);
        builder.finish()
    }

    /// Write these metrics to the HDU `AO_FLAG_METRICS` in an open FITS file:
    /// an image of [`AOFlagMetrics::occupancy_tf`] with FITS axes (frequency,
    /// time).
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the FITS file fails.
    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> FitsResult<()> {
        let (num_times, num_freqs) = self.occupancy_tf.dim();
        let hdu = write_image_hdu(
            fptr,
            "AO_FLAG_METRICS",
            ImageType::Double,
            &[num_times, num_freqs],
            self.occupancy_tf.iter().copied(),
        )?;
        TimeFreqAxes {
            start_time_gps_s: self.start_time_gps_s,
            integration_time_s: self.integration_time_s,
            start_freq_hz: self.start_freq_hz,
            channel_width_hz: self.channel_width_hz,
        }
        .write(fptr, &hdu)
    }
}

/// Cross-correlation delay metrics for short baselines.
///
/// For each cross-correlation baseline no longer than a cutoff, the
/// visibilities are averaged coherently over the unflagged timesteps, and the
/// amplitude of that spectrum is delay transformed between 100 and 3000 ns
/// (see [`mod@crate::delay_transform`]). The metric for each tile is the mean of
/// the delay spectra of the short baselines it is part of.
///
/// The polarisations are in the order of a [`Jones`] matrix: `XX, XY, YX, YY`.
/// The antenna axis is in order of increasing antenna index.
pub struct CrossMetrics {
    /// Mean delay spectrum of the short baselines of each antenna
    /// `(antenna, delay, polarisation)`. `NaN` for an antenna with no short
    /// baselines.
    pub short_ant_delay_pol_adp: Array3<f32>,
    /// Tile names of the selected antennas
    pub antenna_names: Vec<String>,
    /// Antenna indices of the selected antennas, which is the order of the
    /// antennas in the metafits
    pub antenna_ids: Vec<u32>,
    /// Tile IDs of the selected antennas
    pub antenna_nums: Vec<u32>,
    /// Centre frequency of the first channel \[Hz\]
    pub start_freq_hz: f64,
    /// Fine channel width \[Hz\]
    pub channel_width_hz: f64,
    /// GPS time of the first timestep \[s\]
    pub start_time_gps_s: f64,
    /// Integration time of each timestep \[s\]
    pub integration_time_s: f64,
    /// Delay of the first bin of `short_ant_delay_pol_adp` \[ns\]
    pub start_delay_ns: f64,
    /// Width of the bins of `short_ant_delay_pol_adp` \[ns\]
    pub delay_res_ns: f64,
}

/// Computes [`CrossMetrics`] from visibilities which arrive in chunks of
/// timesteps, so that the visibilities of the whole observation don't need to
/// be in memory at once.
///
/// The result is the same as [`CrossMetrics::new_from_metadata`] on all of the
/// visibilities, however they are chunked.
pub struct CrossMetricsBuilder {
    metadata: MetricsContext,
    timestep_flags: Vec<bool>,
    chan_flags: Vec<bool>,
    /// The baseline index of each short baseline, and the positions of its
    /// antennas in `metadata.antennas`.
    short_baselines: Vec<(usize, usize, usize)>,
    /// Sum over time of the finite visibilities of each short baseline
    /// `(short baseline, frequency, polarisation)`
    sum_vis_bfp: Array3<Complex<f64>>,
    /// How many visibilities are in each sum
    count_bfp: Array3<u32>,
    num_timesteps_pushed: usize,
}

impl CrossMetricsBuilder {
    /// Start computing [`CrossMetrics`] for all of the timesteps, channels
    /// and baselines in `metadata`, using baselines up to `baseline_cutoff_m`
    /// metres long.
    pub fn new(
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
        baseline_cutoff_m: f32,
    ) -> Self {
        // metadata.antennas is in the order of the sorted antenna indices.
        let sel_ants = selected_antennas(&metadata.antenna_pairs);
        let ant_pos = |ant_idx: usize| {
            sel_ants
                .binary_search(&ant_idx)
                .expect("antenna is in a pair")
        };
        let cutoff_sq = f64::from(baseline_cutoff_m * baseline_cutoff_m);

        let short_baselines: Vec<_> = metadata
            .antenna_pairs
            .iter()
            .enumerate()
            .filter(|&(_, &(ant1, ant2))| ant1 != ant2)
            .map(|(bl_idx, &(ant1, ant2))| (bl_idx, ant_pos(ant1), ant_pos(ant2)))
            .filter(|&(_, pos1, pos2)| {
                let (meta1, meta2) = (&metadata.antennas[pos1], &metadata.antennas[pos2]);
                let dx = meta1.north_m - meta2.north_m;
                let dy = meta1.east_m - meta2.east_m;
                let dz = meta1.height_m - meta2.height_m;
                dx * dx + dy * dy + dz * dz <= cutoff_sq
            })
            .collect();

        let sum_dim = (short_baselines.len(), chan_flags.len(), 4);
        Self {
            metadata: metadata.clone(),
            timestep_flags: timestep_flags.to_vec(),
            chan_flags: chan_flags.to_vec(),
            short_baselines,
            sum_vis_bfp: Array3::zeros(sum_dim),
            count_bfp: Array3::zeros(sum_dim),
            num_timesteps_pushed: 0,
        }
    }

    /// Add the visibilities `[timestep][channel][baseline]` of the next
    /// timesteps.
    ///
    /// # Panics
    ///
    /// Panics if the channels or baselines don't match the metadata, or if
    /// there are more timesteps than are left.
    pub fn push(&mut self, jones_array_tfb: ArrayView3<Jones<f32>>) {
        let (num_timesteps, num_freqs, num_baselines) = jones_array_tfb.dim();
        let first_timestep = self.num_timesteps_pushed;
        assert!(
            first_timestep + num_timesteps <= self.timestep_flags.len(),
            "too many timesteps"
        );
        assert_eq!(num_freqs, self.chan_flags.len());
        assert_eq!(num_baselines, self.metadata.antenna_pairs.len());

        let timestep_flags = &self.timestep_flags[first_timestep..first_timestep + num_timesteps];
        let chan_flags = &self.chan_flags;
        self.sum_vis_bfp
            .outer_iter_mut()
            .into_par_iter()
            .zip(self.count_bfp.outer_iter_mut())
            .zip(&self.short_baselines)
            .for_each(|((mut sum_fp, mut count_fp), &(bl_idx, _, _))| {
                for t in (0..num_timesteps).filter(|&t| !timestep_flags[t]) {
                    for f in (0..num_freqs).filter(|&f| !chan_flags[f]) {
                        for p in 0..4 {
                            let val = jones_array_tfb[[t, f, bl_idx]][p];
                            if val.re.is_finite() && val.im.is_finite() {
                                sum_fp[[f, p]] +=
                                    Complex::new(f64::from(val.re), f64::from(val.im));
                                count_fp[[f, p]] += 1;
                            }
                        }
                    }
                }
            });
        self.num_timesteps_pushed += num_timesteps;
    }

    /// Compute the metrics.
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - not all of the timesteps have been pushed,
    /// - there are fewer than 2 channels.
    pub fn finish(self) -> CrossMetrics {
        let Self {
            metadata,
            timestep_flags,
            chan_flags,
            short_baselines,
            sum_vis_bfp,
            count_bfp,
            num_timesteps_pushed,
        } = self;
        assert_eq!(
            num_timesteps_pushed,
            timestep_flags.len(),
            "not all timesteps were pushed"
        );
        let num_freqs = chan_flags.len();
        let num_sel_ants = metadata.antennas.len();

        let freqs_arr = Array1::from(metadata.fine_chan_freqs_hz.clone());
        let delay_info = calculate_delay_channels(num_freqs, &freqs_arr, &DELAY_TRANSFORM_CONFIG);
        let num_delays = delay_info.n_delay_channels;

        // The delay spectrum of the amplitude of the time-averaged
        // visibilities of each polarisation of each short baseline, or None if
        // it has no visibilities.
        let delay_spectra: Vec<[Option<Array1<f32>>; 4]> = sum_vis_bfp
            .outer_iter()
            .into_par_iter()
            .zip(count_bfp.outer_iter())
            .map(|(sum_fp, count_fp)| {
                std::array::from_fn(|p| {
                    let spectrum: Array1<f64> = sum_fp
                        .column(p)
                        .iter()
                        .zip(count_fp.column(p))
                        .map(|(sum, &count)| {
                            if count > 0 {
                                (sum / f64::from(count)).norm()
                            } else {
                                0.0
                            }
                        })
                        .collect();
                    if !spectrum.iter().any(|&x| x > 0.0) {
                        return None;
                    }
                    let delay_result = delay_transform(
                        &spectrum.insert_axis(Axis(0)),
                        &freqs_arr,
                        &DELAY_TRANSFORM_CONFIG,
                    )
                    .ok()?;
                    Some(delay_result.delay_spectrum.row(0).mapv(|x| x as f32))
                })
            })
            .collect();

        // The mean of the delay spectra of the short baselines of each antenna
        let mut short_ant_delay_pol_adp = Array3::<f32>::zeros((num_sel_ants, num_delays, 4));
        let mut ant_baseline_counts = Array2::<u32>::zeros((num_sel_ants, 4));
        for (&(_, pos1, pos2), delay_spectra_p) in short_baselines.iter().zip(&delay_spectra) {
            for (p, delay_spectrum) in delay_spectra_p.iter().enumerate() {
                let Some(delay_spectrum) = delay_spectrum else {
                    continue;
                };
                for pos in [pos1, pos2] {
                    let mut acc = short_ant_delay_pol_adp.slice_mut(s![pos, .., p]);
                    acc += delay_spectrum;
                    ant_baseline_counts[[pos, p]] += 1;
                }
            }
        }
        for ((a, p), &count) in ant_baseline_counts.indexed_iter() {
            let mut delay_spectrum = short_ant_delay_pol_adp.slice_mut(s![a, .., p]);
            if count > 0 {
                delay_spectrum /= count as f32;
            } else {
                delay_spectrum.fill(f32::NAN);
            }
        }

        CrossMetrics {
            short_ant_delay_pol_adp,
            antenna_names: metadata
                .antennas
                .iter()
                .map(|ant| ant.tile_name.clone())
                .collect(),
            antenna_ids: metadata.antennas.iter().map(|ant| ant.ant_id).collect(),
            antenna_nums: metadata.antennas.iter().map(|ant| ant.tile_id).collect(),
            start_freq_hz: metadata.fine_chan_freqs_hz[0],
            channel_width_hz: metadata.channel_width_hz(),
            start_time_gps_s: metadata.timestamps_s[0],
            integration_time_s: metadata.integration_time_s(),
            start_delay_ns: delay_info.start_delay_ns,
            delay_res_ns: delay_info.actual_delay_res_ns,
        }
    }
}

impl CrossMetrics {
    /// Create [`CrossMetrics`] for visibilities that were read with the
    /// selection `chunk_vis_sel` of the observation `corr_ctx`, using
    /// baselines up to `baseline_cutoff_m` metres long.
    ///
    /// # Panics
    ///
    /// See [`CrossMetrics::new_from_metadata`].
    pub fn new(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        corr_ctx: &CorrelatorContext,
        chunk_vis_sel: &VisSelection,
        flag_ctx: &FlagContext,
        baseline_cutoff_m: f32,
    ) -> Self {
        let (metadata, timestep_flags, chan_flags) =
            selection_metadata(corr_ctx, chunk_vis_sel, flag_ctx);
        Self::new_from_metadata(
            jones_array_tfb,
            &metadata,
            &timestep_flags,
            &chan_flags,
            baseline_cutoff_m,
        )
    }

    /// Create [`CrossMetrics`] from visibilities `[timestep][channel][baseline]`
    /// and their metadata, using baselines up to `baseline_cutoff_m` metres
    /// long.
    ///
    /// To compute the metrics from visibilities that don't fit in memory, see
    /// [`CrossMetricsBuilder`].
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - there are fewer than 2 channels,
    /// - the dimensions of the visibilities or flags don't match `metadata`.
    pub fn new_from_metadata(
        jones_array_tfb: ArrayView3<Jones<f32>>,
        metadata: &MetricsContext,
        timestep_flags: &[bool],
        chan_flags: &[bool],
        baseline_cutoff_m: f32,
    ) -> Self {
        let (num_timesteps, num_freqs, _) = jones_array_tfb.dim();
        let mut builder = CrossMetricsBuilder::new(
            metadata,
            &timestep_flags[..num_timesteps],
            &chan_flags[..num_freqs],
            baseline_cutoff_m,
        );
        builder.push(jones_array_tfb);
        builder.finish()
    }

    /// Write these metrics to HDUs in an open FITS file: `CROSS_DELAY_POL=<pol>`
    /// for each polarisation, which is an image of
    /// [`CrossMetrics::short_ant_delay_pol_adp`] with FITS axes (delay,
    /// antenna).
    ///
    /// The HDUs are in the order of the FITS Stokes axis: `XX, YY, XY, YX`.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the FITS file fails.
    pub fn save_to_fits(&self, fptr: &mut FitsFile) -> FitsResult<()> {
        let (num_ants, num_delays, _) = self.short_ant_delay_pol_adp.dim();

        for pol_name in FITS_POL_NAMES {
            // the array is in the order of a Jones matrix
            let p = JONES_POL_NAMES
                .iter()
                .position(|&name| name == pol_name)
                .expect("all pols are in a Jones matrix");
            let hdu = write_image_hdu(
                fptr,
                &format!("CROSS_DELAY_POL={pol_name}"),
                ImageType::Float,
                &[num_ants, num_delays],
                self.short_ant_delay_pol_adp
                    .slice(s![.., .., p])
                    .iter()
                    .copied(),
            )?;
            FitsAxis {
                ctype: "DELAY",
                crval: self.start_delay_ns,
                cdelt: None,
                cunit: Some("ns"),
            }
            .write(fptr, &hdu, 1)?;
            hdu.write_key(fptr, "DELAY_RES", self.delay_res_ns)?;
            hdu.write_key(fptr, "CTYPE2", "ANTENNA")?;
            hdu.write_key(fptr, "POL", pol_name)?;
            hdu.write_key(fptr, "N_ANTS", num_ants as u32)?;
            hdu.write_key(fptr, "TELESCOP", "MWA")?;
            hdu.write_key(fptr, "INSTRUME", "CROSS_METRICS")?;
            hdu.write_key(fptr, "ORIGIN", "Birli")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod test_common {
    use super::{AntennaMetadata, MetricsContext};
    use crate::marlu::{
        fitsio::FitsFile,
        ndarray::{Array3, ArrayD, IxDyn},
        num_complex::Complex,
        Jones,
    };
    use std::f64::consts::TAU;
    use tempfile::NamedTempFile;

    /// Metadata for `num_times` timesteps of 2 s from GPS time 1e9,
    /// `num_freqs` channels of 40 kHz from 150 MHz, and all the baselines
    /// (with autos) of antennas with indices `ant_idxs`, 10 m apart on a
    /// north-south line in order of index.
    pub(super) fn metrics_context(
        num_times: usize,
        num_freqs: usize,
        ant_idxs: &[usize],
    ) -> MetricsContext {
        let mut antenna_pairs = vec![];
        for (i, &ant1) in ant_idxs.iter().enumerate() {
            for &ant2 in &ant_idxs[i..] {
                antenna_pairs.push((ant1, ant2));
            }
        }
        let mut sorted_ant_idxs = ant_idxs.to_vec();
        sorted_ant_idxs.sort_unstable();
        MetricsContext {
            antennas: sorted_ant_idxs
                .iter()
                .map(|&idx| AntennaMetadata {
                    tile_name: format!("Tile{idx:03}"),
                    tile_id: 1000 + idx as u32,
                    ant_id: idx as u32,
                    east_m: 0.5,
                    north_m: 10.0 * idx as f64,
                    height_m: 1.5,
                    rec_number: 1 + idx as u32 / 8,
                    rec_slot_number: 1 + idx as u32 % 8,
                    rec_type: "RRI".into(),
                    cable_flavour: format!("RG6_{}", 90 + 60 * (idx % 2)),
                    has_whitening_filter: idx % 2 == 0,
                })
                .collect(),
            fine_chan_freqs_hz: (0..num_freqs).map(|f| 150e6 + 40e3 * f as f64).collect(),
            timestamps_s: (0..num_times).map(|t| 1e9 + 2.0 * t as f64).collect(),
            antenna_pairs,
        }
    }

    /// Visibilities for `ctx` which are different for every timestep, channel,
    /// baseline and polarisation.
    pub(super) fn synthetic_vis(ctx: &MetricsContext) -> Array3<Jones<f32>> {
        let dims = (
            ctx.timestamps_s.len(),
            ctx.fine_chan_freqs_hz.len(),
            ctx.antenna_pairs.len(),
        );
        Array3::from_shape_fn(dims, |(t, f, b)| {
            let (ant1, ant2) = ctx.antenna_pairs[b];
            Jones::from([0, 1, 2, 3].map(|p| {
                // a bandpass with a ripple whose period depends on the baseline
                // and polarisation, and which varies in time
                let ripple = (TAU * f as f64 / (4.0 + (ant1 + 2 * ant2 + p) as f64)).cos();
                let amp = (1 + ant1 + ant2) as f64 * (2.0 + ripple) + 0.1 * ((t * 7) % 5) as f64;
                let phase = 0.3 * (p + b) as f64 + 0.01 * (t * f) as f64;
                Complex::from_polar(amp as f32, phase as f32)
            }))
        })
    }

    /// Write something to a temporary FITS file with `save`, then open it.
    pub(super) fn save_and_open(
        save: impl FnOnce(&mut FitsFile) -> super::FitsResult<()>,
    ) -> (NamedTempFile, FitsFile) {
        let tmp_file = NamedTempFile::new().unwrap();
        std::fs::remove_file(tmp_file.path()).unwrap();
        {
            let mut fptr = FitsFile::create(tmp_file.path()).open().unwrap();
            save(&mut fptr).unwrap();
        }
        let fptr = FitsFile::open(tmp_file.path()).unwrap();
        (tmp_file, fptr)
    }

    /// The names of the HDUs in a FITS file, after the primary HDU.
    pub(super) fn hdu_names(fptr: &mut FitsFile) -> Vec<String> {
        let num_hdus = fptr.iter().count();
        (1..num_hdus)
            .map(|i| {
                let hdu = fptr.hdu(i).unwrap();
                hdu.read_key::<String>(fptr, "EXTNAME").unwrap()
            })
            .collect()
    }

    /// Read an image HDU into an array with the dimensions it was written with.
    pub(super) fn read_image(fptr: &mut FitsFile, extname: &str) -> ArrayD<f64> {
        let hdu = fptr.hdu(extname).unwrap();
        // the shape is in the order of the array: the reverse of the FITS axes
        let num_axes: i64 = hdu.read_key(fptr, "NAXIS").unwrap();
        let shape: Vec<usize> = (1..=num_axes)
            .rev()
            .map(|axis| hdu.read_key::<i64>(fptr, &format!("NAXIS{axis}")).unwrap() as usize)
            .collect();
        let data: Vec<f64> = hdu.read_image(fptr).unwrap();
        ArrayD::from_shape_vec(IxDyn(&shape), data).unwrap()
    }

    /// Assert that the string keys of an HDU are as expected.
    pub(super) fn assert_str_keys(fptr: &mut FitsFile, extname: &str, expected: &[(&str, &str)]) {
        let hdu = fptr.hdu(extname).unwrap();
        for &(key, value) in expected {
            let observed: String = hdu.read_key(fptr, key).unwrap();
            assert_eq!(observed, value, "{extname} {key}");
        }
    }

    /// Assert that the numeric keys of an HDU are as expected.
    pub(super) fn assert_f64_keys(fptr: &mut FitsFile, extname: &str, expected: &[(&str, f64)]) {
        let hdu = fptr.hdu(extname).unwrap();
        for &(key, value) in expected {
            let observed: f64 = hdu.read_key(fptr, key).unwrap();
            assert!(
                (observed - value).abs() <= 1e-9 * value.abs(),
                "{extname} {key}: {observed} != {value}"
            );
        }
    }

    /// Assert that two arrays are the same to within a relative tolerance,
    /// with NaNs in the same places.
    pub(super) fn assert_close<'a>(
        observed: impl IntoIterator<Item = &'a f32>,
        expected: impl IntoIterator<Item = f64>,
        rel_tol: f64,
    ) {
        let observed: Vec<f64> = observed.into_iter().map(|&x| f64::from(x)).collect();
        let expected: Vec<f64> = expected.into_iter().collect();
        assert_eq!(observed.len(), expected.len());
        for (i, (obs, exp)) in observed.iter().zip(&expected).enumerate() {
            assert!(
                (obs.is_nan() && exp.is_nan()) || (obs - exp).abs() <= rel_tol * exp.abs().max(1.0),
                "element {i}: {obs} != {exp}"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod context_tests {
    use super::{hyperdrive_to_fits_stokes, selected_antennas, MetricsContext};
    use crate::marlu::{mwalib::CorrelatorContext, VisSelection};

    fn get_mwax_context() -> CorrelatorContext {
        CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &[
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_hyperdrive_to_fits_stokes() {
        // XX XY YX YY -> XX YY XY YX
        let jones_pols = super::JONES_POL_NAMES;
        let fits_pols = super::FITS_POL_NAMES;
        for (p, jones_pol) in jones_pols.into_iter().enumerate() {
            assert_eq!(fits_pols[hyperdrive_to_fits_stokes(p)], jones_pol);
        }
    }

    #[test]
    #[should_panic(expected = "Invalid pol: 4")]
    fn test_hyperdrive_to_fits_stokes_invalid() {
        hyperdrive_to_fits_stokes(4);
    }

    #[test]
    fn test_selected_antennas() {
        assert_eq!(
            selected_antennas(&[(5, 5), (5, 2), (2, 2), (7, 2)]),
            vec![2, 5, 7]
        );
        assert!(selected_antennas(&[]).is_empty());
    }

    #[test]
    fn test_from_mwalib() {
        let corr_ctx = get_mwax_context();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.timestep_range = 1..3;
        vis_sel.coarse_chan_range = 1..2;

        let ctx = MetricsContext::from_mwalib(&corr_ctx, &vis_sel);

        // every antenna is selected, in order
        let meta_ctx = &corr_ctx.metafits_context;
        assert_eq!(ctx.antennas.len(), meta_ctx.num_ants);
        assert_eq!(
            ctx.antenna_pairs,
            vis_sel.get_ant_pairs(&corr_ctx.metafits_context)
        );
        for (ant, meta_ant) in ctx.antennas.iter().zip(&meta_ctx.antennas) {
            assert_eq!(ant.tile_name, meta_ant.tile_name);
            assert_eq!(ant.tile_id, meta_ant.tile_id);
            assert_eq!(ant.ant_id, meta_ant.ant);
            assert_eq!(ant.north_m.to_bits(), meta_ant.north_m.to_bits());
            assert_eq!(ant.east_m.to_bits(), meta_ant.east_m.to_bits());
            assert_eq!(ant.height_m.to_bits(), meta_ant.height_m.to_bits());
            assert_eq!(ant.rec_number, meta_ant.rfinput_x.rec_number);
            assert_eq!(ant.rec_slot_number, meta_ant.rfinput_x.rec_slot_number);
            assert_eq!(ant.rec_type, meta_ant.rfinput_x.rec_type.to_string());
            assert_eq!(ant.cable_flavour, meta_ant.rfinput_x.flavour);
            assert_eq!(
                ant.has_whitening_filter,
                meta_ant.rfinput_x.has_whitening_filter
            );
        }

        // the selected timesteps and coarse channel
        assert_eq!(
            ctx.timestamps_s,
            corr_ctx.timesteps[1..3]
                .iter()
                .map(|ts| ts.gps_time_ms as f64 / 1000.0)
                .collect::<Vec<_>>()
        );
        let fine_chans_per_coarse = meta_ctx.num_corr_fine_chans_per_coarse;
        assert_eq!(ctx.fine_chan_freqs_hz.len(), fine_chans_per_coarse);
        assert_eq!(
            ctx.fine_chan_freqs_hz,
            meta_ctx.metafits_fine_chan_freqs_hz[fine_chans_per_coarse..2 * fine_chans_per_coarse]
        );

        let int_time_s = meta_ctx.corr_int_time_ms as f64 / 1000.0;
        assert!((ctx.integration_time_s() - int_time_s).abs() < 1e-9);
        let fine_chan_width_hz = f64::from(meta_ctx.corr_fine_chan_width_hz);
        assert!((ctx.channel_width_hz() - fine_chan_width_hz).abs() < 1e-9);
    }

    /// The antennas are in order of index, whatever the order of the baselines.
    #[test]
    fn test_from_mwalib_antenna_selection() {
        let corr_ctx = get_mwax_context();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.retain_antennas(&corr_ctx.metafits_context, &[1]);
        vis_sel.timestep_range = 0..1;

        let ctx = MetricsContext::from_mwalib(&corr_ctx, &vis_sel);

        assert_eq!(ctx.antenna_pairs, vec![(1, 1)]);
        assert_eq!(ctx.antennas.len(), 1);
        assert_eq!(ctx.antennas[0].ant_id, 1);
        assert_eq!(
            ctx.antennas[0].tile_name,
            corr_ctx.metafits_context.antennas[1].tile_name
        );
        // with one timestep, the integration time can't be measured
        assert!((ctx.integration_time_s() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_single_channel_width() {
        let mut ctx = super::test_common::metrics_context(2, 1, &[0]);
        assert!((ctx.channel_width_hz() - 1.0).abs() < 1e-9);
        ctx.fine_chan_freqs_hz.push(150.04e6);
        assert!((ctx.channel_width_hz() - 40e3).abs() < 1e-3);
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod autometrics_tests {
    use super::{
        test_common::{
            assert_close, assert_f64_keys, assert_str_keys, hdu_names, metrics_context, read_image,
            save_and_open, synthetic_vis,
        },
        AutoMetrics, FITS_POL_NAMES,
    };
    use crate::{
        delay_transform::{delay_transform, DelayTransformConfig},
        marlu::{
            mwalib::CorrelatorContext,
            ndarray::{s, Array1, Array2, Array3},
            Jones, VisSelection,
        },
        FlagContext,
    };

    /// index into a Jones matrix of each polarisation in FITS order
    const FITS_TO_JONES: [usize; 4] = [0, 3, 1, 2];

    /// The values are right, with antennas that are not in order, and a
    /// flagged timestep and channel.
    #[test]
    fn test_autometrics_values() {
        let (num_times, num_freqs) = (5, 16);
        // antenna 7 comes before antenna 2 in the baselines
        let ctx = metrics_context(num_times, num_freqs, &[7, 2, 4]);
        let vis = synthetic_vis(&ctx);
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[3] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[5] = true;

        let auto = AutoMetrics::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags);

        // the antennas are in order of index
        assert_eq!(auto.antenna_names, ["Tile002", "Tile004", "Tile007"]);
        assert_eq!(auto.antenna_ids, [2, 4, 7]);
        assert_eq!(auto.antenna_nums, [1002, 1004, 1007]);
        assert_eq!(auto.rx_numbers, [1, 1, 1]);
        assert_eq!(auto.rx_slots, [3, 5, 8]);
        assert_eq!(auto.rx_types, ["RRI", "RRI", "RRI"]);
        assert_eq!(auto.cable_flavours, ["RG6_90", "RG6_90", "RG6_150"]);
        assert_eq!(auto.whitening_filters, [true, true, false]);
        for (position, north_m) in auto.antenna_positions.iter().zip([20.0, 40.0, 70.0]) {
            assert_eq!((position.n, position.e, position.h), (north_m, 0.5, 1.5));
        }
        assert_eq!(auto.start_freq_hz.to_bits(), 150e6_f64.to_bits());
        assert_eq!(auto.channel_width_hz.to_bits(), 40e3_f64.to_bits());
        assert_eq!(auto.start_time_gps_s.to_bits(), 1e9_f64.to_bits());
        assert_eq!(auto.integration_time_s.to_bits(), 2.0_f64.to_bits());

        assert_eq!(auto.auto_sub_aptf.dim(), (3, 4, num_times, num_freqs));
        assert_eq!(auto.auto_spectrum_afp.dim(), (3, num_freqs, 4));
        assert_eq!(auto.auto_coeffs_apo.dim(), (3, 4, 3));

        let unflagged_times = [0, 1, 2, 4];
        for (a, ant_idx) in [2, 4, 7].into_iter().enumerate() {
            let bl_idx = ctx
                .antenna_pairs
                .iter()
                .position(|&pair| pair == (ant_idx, ant_idx))
                .unwrap();
            for (fits_p, jones_p) in FITS_TO_JONES.into_iter().enumerate() {
                let amp = |t: usize, f: usize| f64::from(vis[[t, f, bl_idx]][jones_p].norm());
                for f in 0..num_freqs {
                    let spectrum = auto.auto_spectrum_afp[[a, f, fits_p]];
                    let sub = auto.auto_sub_aptf.slice(s![a, fits_p, .., f]);
                    if chan_flags[f] {
                        assert!(spectrum.is_nan());
                        assert!(sub.iter().all(|x| x.is_nan()));
                        continue;
                    }
                    // the bandpass is the mean over unflagged timesteps
                    let mean = unflagged_times.iter().map(|&t| amp(t, f)).sum::<f64>() / 4.0;
                    assert_close([&spectrum], [mean], 1e-6);
                    // and it is subtracted from the amplitudes
                    let expected = (0..num_times).map(|t| {
                        if timestep_flags[t] {
                            f64::NAN
                        } else {
                            amp(t, f) - mean
                        }
                    });
                    assert_close(sub.iter(), expected, 1e-5);
                }
            }
        }

        // the polynomial is a fit to the bandpass: its residuals are no
        // bigger than the ripple.
        for ((a, f, p), &spectrum) in auto.auto_spectrum_afp.indexed_iter() {
            if chan_flags[f] {
                continue;
            }
            let freq_hz = ctx.fine_chan_freqs_hz[f];
            let fit: f64 = (0..3)
                .map(|k| f64::from(auto.auto_coeffs_apo[[a, p, k]]) * freq_hz.powi(k as i32))
                .sum();
            let ant_idx = [2, 4, 7][a] as f64;
            assert!(
                (fit - f64::from(spectrum)).abs() < 1.5 * (1.0 + 2.0 * ant_idx),
                "a={a} f={f} p={p} fit={fit} spectrum={spectrum}"
            );
        }

        // the delay spectrum is the delay transform of the bandpass
        let config = DelayTransformConfig {
            min_delay_ns: 100.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1.0,
        };
        let freqs_hz = Array1::from(ctx.fine_chan_freqs_hz.clone());
        for p in 0..4 {
            let spectrum_af: Array2<f64> =
                auto.auto_spectrum_afp.slice(s![.., .., p]).mapv(f64::from);
            let expected = delay_transform(&spectrum_af, &freqs_hz, &config).unwrap();
            assert_close(
                auto.auto_delay_afp.slice(s![.., .., p]).iter(),
                expected.delay_spectrum.iter().copied(),
                1e-6,
            );
            assert!((auto.start_delay_ns - expected.delays_ns[0]).abs() < 1e-9);
            assert!(
                (auto.delay_res_ns - (expected.delays_ns[1] - expected.delays_ns[0])).abs() < 1e-9
            );
        }
        assert!(auto.start_delay_ns >= 100.0 && auto.start_delay_ns < 101.0);
        assert!(auto.delay_res_ns <= 1.0 && auto.delay_res_ns > 0.9);
    }

    /// Only the antennas with autocorrelations have metrics, and they get the
    /// right metadata.
    #[test]
    fn test_autometrics_without_some_autos() {
        let mut ctx = metrics_context(2, 4, &[1, 3, 6]);
        let vis = synthetic_vis(&ctx);
        // leave out the autocorrelation of antenna 3
        let keep: Vec<usize> = (0..ctx.antenna_pairs.len())
            .filter(|&b| ctx.antenna_pairs[b] != (3, 3))
            .collect();
        let vis = vis.select(crate::marlu::ndarray::Axis(2), &keep);
        ctx.antenna_pairs = keep.iter().map(|&b| ctx.antenna_pairs[b]).collect();

        let auto = AutoMetrics::new_from_metadata(vis.view(), &ctx, &[false; 2], &[false; 4]);

        assert_eq!(auto.antenna_names, ["Tile001", "Tile006"]);
        assert_eq!(auto.antenna_ids, [1, 6]);
        assert_eq!(auto.auto_spectrum_afp.dim(), (2, 4, 4));
        // antenna 6's auto is the last baseline
        let last_bl = ctx.antenna_pairs.len() - 1;
        let amp = (vis[[0, 2, last_bl]][0].norm() + vis[[1, 2, last_bl]][0].norm()) / 2.0;
        assert_close([&auto.auto_spectrum_afp[[1, 2, 0]]], [f64::from(amp)], 1e-6);
    }

    #[test]
    #[should_panic(expected = "all timesteps are flagged")]
    fn test_autometrics_all_timesteps_flagged() {
        let ctx = metrics_context(2, 4, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        AutoMetrics::new_from_metadata(vis.view(), &ctx, &[true; 2], &[false; 4]);
    }

    #[test]
    #[should_panic(expected = "all channels are flagged")]
    fn test_autometrics_all_channels_flagged() {
        let ctx = metrics_context(2, 4, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        AutoMetrics::new_from_metadata(vis.view(), &ctx, &[false; 2], &[true; 4]);
    }

    /// The polynomial coefficients are NaN when there aren't enough unflagged
    /// channels to fit.
    #[test]
    fn test_autometrics_too_few_channels_to_fit() {
        let ctx = metrics_context(2, 4, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        let chan_flags = [false, true, true, false];
        let auto = AutoMetrics::new_from_metadata(vis.view(), &ctx, &[false; 2], &chan_flags);
        assert!(auto.auto_coeffs_apo.iter().all(|c| c.is_nan()));
        assert!(auto.auto_delay_afp.iter().all(|d| d.is_finite()));
    }

    /// `AutoMetrics::new` uses the antennas, flags, times and frequencies of
    /// the selection.
    #[test]
    fn test_autometrics_new_with_selection() {
        let corr_ctx = CorrelatorContext::new(
            "tests/data/1119683928_picket/1119683928.metafits",
            &["tests/data/1119683928_picket/1119683928_20150630071834_gpubox01_00.fits"],
        )
        .unwrap();
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let meta_ants = &corr_ctx.metafits_context.antennas;

        // antennas which are not at the start, sequential or in order.
        for ant_idxs in [vec![1, 2], vec![0, 1], vec![5], vec![20, 4, 19]] {
            let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
            vis_sel.retain_antennas(&corr_ctx.metafits_context, &ant_idxs);
            let mut sorted_ant_idxs = ant_idxs.clone();
            sorted_ant_idxs.sort_unstable();
            let num_ants = ant_idxs.len();

            let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
            flag_ctx.fine_chan_flags[3] = true;
            let (num_times, num_freqs, num_baselines) = vis_sel.get_shape(fine_chans_per_coarse);
            let jones_array =
                Array3::from_shape_fn((num_times, num_freqs, num_baselines), |(_, _, b)| {
                    Jones::identity() * (1.0 + b as f32)
                });

            let auto = AutoMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

            assert_eq!(
                auto.auto_sub_aptf.dim(),
                (num_ants, 4, num_times, num_freqs)
            );
            assert_eq!(auto.auto_spectrum_afp.dim(), (num_ants, num_freqs, 4));
            assert_eq!(auto.auto_coeffs_apo.dim(), (num_ants, 4, 3));
            assert_eq!(auto.auto_delay_afp.dim().0, num_ants);
            for (a, &ant_idx) in sorted_ant_idxs.iter().enumerate() {
                assert_eq!(auto.antenna_names[a], meta_ants[ant_idx].tile_name);
                assert_eq!(auto.antenna_ids[a], ant_idx as u32);
                assert_eq!(auto.antenna_nums[a], meta_ants[ant_idx].tile_id);
                // XX of the autocorrelation of this antenna
                let bl_idx = vis_sel
                    .get_ant_pairs(&corr_ctx.metafits_context)
                    .iter()
                    .position(|&pair| pair == (ant_idx, ant_idx))
                    .unwrap();
                assert_close(
                    [&auto.auto_spectrum_afp[[a, 0, 0]]],
                    [1.0 + bl_idx as f64],
                    1e-6,
                );
                // the flagged fine channel
                assert!(auto.auto_spectrum_afp[[a, 3, 0]].is_nan());
            }
            assert_eq!(
                auto.start_freq_hz.to_bits(),
                corr_ctx.metafits_context.metafits_fine_chan_freqs_hz[0].to_bits()
            );
            assert_eq!(
                auto.start_time_gps_s.to_bits(),
                (corr_ctx.timesteps[vis_sel.timestep_range.start].gps_time_ms as f64 / 1000.0)
                    .to_bits()
            );
        }
    }

    #[test]
    fn test_autometrics_save_to_fits() {
        let (num_times, num_freqs) = (3, 8);
        let ctx = metrics_context(num_times, num_freqs, &[4, 1]);
        let vis = synthetic_vis(&ctx);
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[0] = true;
        let auto = AutoMetrics::new_from_metadata(vis.view(), &ctx, &[false; 3], &chan_flags);
        let num_delays = auto.auto_delay_afp.dim().1;

        let (_tmp_file, mut fptr) = save_and_open(|fptr| auto.save_to_fits(fptr));

        assert_eq!(
            hdu_names(&mut fptr),
            [
                "AUTO_SUB_ANT=Tile001",
                "AUTO_SUB_ANT=Tile004",
                "AUTO_POL=XX",
                "AUTO_POL=YY",
                "AUTO_POL=XY",
                "AUTO_POL=YX",
                "AUTO_COEFFS",
                "AUTO_DELAY_POL=XX",
                "AUTO_DELAY_POL=YY",
                "AUTO_DELAY_POL=XY",
                "AUTO_DELAY_POL=YX",
            ]
        );

        for (a, extname) in ["AUTO_SUB_ANT=Tile001", "AUTO_SUB_ANT=Tile004"]
            .into_iter()
            .enumerate()
        {
            let image = read_image(&mut fptr, extname);
            assert_eq!(image.shape(), [4, num_times, num_freqs]);
            assert_close(
                auto.auto_sub_aptf.slice(s![a, .., .., ..]).iter(),
                image.iter().copied(),
                0.0,
            );
            // the axes are in FITS order: the last dimension is axis 1.
            assert_str_keys(
                &mut fptr,
                extname,
                &[
                    ("CTYPE1", "FREQ"),
                    ("CUNIT1", "Hz"),
                    ("CTYPE2", "TIME"),
                    ("CUNIT2", "s"),
                    ("CTYPE3", "STOKES"),
                    ("ANTNAME", &auto.antenna_names[a]),
                    ("RX_TYPE", "RRI"),
                    ("CABLE_FLAVOUR", &auto.cable_flavours[a]),
                ],
            );
            assert_f64_keys(
                &mut fptr,
                extname,
                &[
                    ("NAXIS1", num_freqs as f64),
                    ("NAXIS2", num_times as f64),
                    ("NAXIS3", 4.0),
                    ("CRVAL1", 150e6),
                    ("CRPIX1", 1.0),
                    ("CHAN_WIDTH", 40e3),
                    ("CRVAL2", 1e9),
                    ("CDELT2", 2.0),
                    ("CRPIX2", 1.0),
                    ("CRVAL3", -5.0),
                    ("CDELT3", -1.0),
                    ("CRPIX3", 1.0),
                    ("ANT_ID", f64::from(auto.antenna_ids[a])),
                    ("ANT_NUM", f64::from(auto.antenna_nums[a])),
                    (
                        "WHITENING_FILTER",
                        f64::from(u8::from(auto.whitening_filters[a])),
                    ),
                    ("RX_NUMBER", f64::from(auto.rx_numbers[a])),
                    ("RX_SLOT", f64::from(auto.rx_slots[a])),
                    ("NORTH", auto.antenna_positions[a].n),
                    ("EAST", auto.antenna_positions[a].e),
                    ("HEIGHT", auto.antenna_positions[a].h),
                ],
            );
        }

        for (p, pol_name) in FITS_POL_NAMES.into_iter().enumerate() {
            let extname = format!("AUTO_POL={pol_name}");
            let image = read_image(&mut fptr, &extname);
            assert_eq!(image.shape(), [2, num_freqs]);
            assert_close(
                auto.auto_spectrum_afp.slice(s![.., .., p]).iter(),
                image.iter().copied(),
                0.0,
            );
            assert_str_keys(
                &mut fptr,
                &extname,
                &[
                    ("CTYPE1", "FREQ"),
                    ("CUNIT1", "Hz"),
                    ("CTYPE2", "BASELINE"),
                    ("POL", pol_name),
                    ("TELESCOP", "MWA"),
                ],
            );
            assert_f64_keys(
                &mut fptr,
                &extname,
                &[("CRVAL1", 150e6), ("CHAN_WIDTH", 40e3), ("N_ANTS", 2.0)],
            );
            // XX of the autocorrelation of antenna 4, which is the first
            // baseline, is the second row.
            if pol_name == "XX" {
                let mean_amp = (0..num_times)
                    .map(|t| f64::from(vis[[t, 2, 0]][0].norm()))
                    .sum::<f64>()
                    / num_times as f64;
                assert!((image[[1, 2]] - mean_amp).abs() < 1e-5);
            }

            let extname = format!("AUTO_DELAY_POL={pol_name}");
            let image = read_image(&mut fptr, &extname);
            assert_eq!(image.shape(), [2, num_delays]);
            assert_close(
                auto.auto_delay_afp.slice(s![.., .., p]).iter(),
                image.iter().copied(),
                0.0,
            );
            assert_str_keys(
                &mut fptr,
                &extname,
                &[
                    ("CTYPE1", "DELAY"),
                    ("CUNIT1", "ns"),
                    ("CTYPE2", "BASELINE"),
                    ("POL", pol_name),
                ],
            );
            assert_f64_keys(
                &mut fptr,
                &extname,
                &[
                    ("CRVAL1", auto.start_delay_ns),
                    ("DELAY_RES", auto.delay_res_ns),
                    ("N_ANTS", 2.0),
                ],
            );
        }

        let hdu = fptr.hdu("AUTO_COEFFS").unwrap();
        let names: Vec<String> = hdu.read_col(&mut fptr, "ANT_NAME").unwrap();
        assert_eq!(names, auto.antenna_names);
        let ant_ids: Vec<u32> = hdu.read_col(&mut fptr, "ANT_ID").unwrap();
        assert_eq!(ant_ids, auto.antenna_ids);
        let ant_nums: Vec<u32> = hdu.read_col(&mut fptr, "ANT_NUM").unwrap();
        assert_eq!(ant_nums, auto.antenna_nums);
        let flavours: Vec<String> = hdu.read_col(&mut fptr, "CABLE_FLAVOUR").unwrap();
        assert_eq!(flavours, auto.cable_flavours);
        let whitening_filters: Vec<i32> = hdu.read_col(&mut fptr, "WHITENING_FILTER").unwrap();
        assert_eq!(whitening_filters, [0, 1]);
        let rx_numbers: Vec<u32> = hdu.read_col(&mut fptr, "RX_NUMBER").unwrap();
        assert_eq!(rx_numbers, auto.rx_numbers);
        let rx_slots: Vec<u32> = hdu.read_col(&mut fptr, "RX_SLOT").unwrap();
        assert_eq!(rx_slots, auto.rx_slots);
        let rx_types: Vec<String> = hdu.read_col(&mut fptr, "RX_TYPE").unwrap();
        assert_eq!(rx_types, auto.rx_types);
        let north_m: Vec<f64> = hdu.read_col(&mut fptr, "NORTH").unwrap();
        assert_eq!(north_m, [10.0, 40.0]);
        let east_m: Vec<f64> = hdu.read_col(&mut fptr, "EAST").unwrap();
        assert_eq!(east_m, [0.5, 0.5]);
        let height_m: Vec<f64> = hdu.read_col(&mut fptr, "HEIGHT").unwrap();
        assert_eq!(height_m, [1.5, 1.5]);
        for (p, pol_name) in FITS_POL_NAMES.into_iter().enumerate() {
            // 3 coefficients for each of the 2 antennas
            let coeffs: Vec<f64> = hdu
                .read_col(&mut fptr, &format!("COEFFS_{pol_name}"))
                .unwrap();
            assert_close(
                auto.auto_coeffs_apo.slice(s![.., p, ..]).iter(),
                coeffs,
                0.0,
            );
        }
        assert_eq!(hdu.read_key::<i64>(&mut fptr, "POLY_ORD").unwrap(), 2);
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod ssins_fits_tests {
    use super::{
        test_common::{
            assert_close, assert_f64_keys, assert_str_keys, hdu_names, metrics_context, read_image,
            save_and_open, synthetic_vis,
        },
        JONES_POL_NAMES, SSINS,
    };
    use crate::{
        marlu::{
            mwalib::CorrelatorContext,
            ndarray::{s, Array3},
            VisSelection,
        },
        FlagContext,
    };

    #[test]
    #[should_panic(expected = "SSINS requires at least 2 timesteps")]
    fn test_ssins_needs_two_timesteps() {
        let ctx = metrics_context(1, 4, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        let flags = Array3::from_elem(vis.dim(), false);
        SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &[false], &[false; 4]);
    }

    /// Flagged timesteps and channels are left out.
    #[test]
    fn test_ssins_timestep_and_chan_flags() {
        let (num_times, num_freqs) = (6, 4);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 2]);
        let vis = synthetic_vis(&ctx);
        let flags = Array3::from_elem(vis.dim(), false);
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[2] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[1] = true;

        let ssins =
            SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &timestep_flags, &chan_flags);

        assert_eq!(ssins.zscore.dim(), (num_times - 1, num_freqs, 4));
        assert_eq!(ssins.num_baselines, 3);
        for ((t, f), &num_baselines) in ssins.num_baselines_tf.indexed_iter() {
            // differences 1 and 2 use the flagged timestep 2
            let flagged = t == 1 || t == 2 || f == 1;
            assert_eq!(num_baselines, if flagged { 0 } else { 3 });
            for p in 0..4 {
                assert_eq!(ssins.zscore[[t, f, p]].is_nan(), flagged);
            }
        }
        assert!(ssins.diff_mean_amp_fp.row(1).iter().all(|x| x.is_nan()));
        assert!(ssins.diff_mean_amp_fp.row(0).iter().all(|x| *x > 0.0));
        assert!(ssins.flag_array.iter().all(|&flag| !flag));
        // the first difference is between the first two timesteps
        assert_eq!(ssins.start_time_gps_s.to_bits(), (1e9_f64 + 1.0).to_bits());
        assert_eq!(ssins.integration_time_s.to_bits(), 2.0_f64.to_bits());
        assert_eq!(ssins.start_freq_hz.to_bits(), 150e6_f64.to_bits());
        assert_eq!(ssins.channel_width_hz.to_bits(), 40e3_f64.to_bits());
    }

    /// `SSINS::new` gives the same thing as `SSINS::new_from_metadata` with
    /// the metadata and flags of the selection.
    #[test]
    fn test_ssins_new_with_selection() {
        let corr_ctx = CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &[
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
            ],
        )
        .unwrap();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        flag_ctx.timestep_flags[3] = true;
        flag_ctx.fine_chan_flags[1] = true;
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        crate::io::read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        let ssins = SSINS::new(
            jones_array.view(),
            flag_array.view(),
            &corr_ctx,
            &vis_sel,
            &flag_ctx,
        );

        let (num_times, num_freqs, _) = jones_array.dim();
        let (ctx, timestep_flags, chan_flags) =
            super::selection_metadata(&corr_ctx, &vis_sel, &flag_ctx);
        assert_eq!(timestep_flags, [false, false, false, true]);
        assert_eq!(chan_flags, [false, true, false, true]);
        let expected = SSINS::new_from_metadata(
            jones_array.view(),
            flag_array.view(),
            &ctx,
            &timestep_flags,
            &chan_flags,
        );
        assert_eq!(ssins.zscore.dim(), (num_times - 1, num_freqs, 4));
        assert_close(
            ssins.zscore.iter(),
            expected.zscore.iter().map(|&z| f64::from(z)),
            0.0,
        );
        assert_eq!(ssins.num_baselines_tf, expected.num_baselines_tf);
        // 2 tiles, so there is one cross-correlation
        assert_eq!(ssins.num_baselines, 1);
        assert!(ssins
            .zscore
            .slice(s![..2, 0, ..])
            .iter()
            .all(|z| z.is_finite()));
        assert!(ssins.zscore.slice(s![2, .., ..]).iter().all(|z| z.is_nan()));
        assert!(ssins.zscore.slice(s![.., 1, ..]).iter().all(|z| z.is_nan()));
    }

    #[test]
    fn test_ssins_save_to_fits() {
        let (num_times, num_freqs) = (4, 6);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 2]);
        let vis = synthetic_vis(&ctx);
        let mut flags = Array3::from_elem(vis.dim(), false);
        // flag a cross-correlation
        flags[[1, 2, 1]] = true;
        let mut ssins =
            SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &[false; 4], &[false; 6]);
        ssins.flag_array[[2, 4]] = true;

        let (_tmp_file, mut fptr) = save_and_open(|fptr| ssins.save_to_fits(fptr));

        assert_eq!(
            hdu_names(&mut fptr),
            [
                "SSINS_POL=XX",
                "SSINS_POL=XY",
                "SSINS_POL=YX",
                "SSINS_POL=YY",
                "SSINS_FLAGS",
                "SSINS_N_BL",
                "SSINS_DIFF_MEAN_AMP_FP",
            ]
        );

        let time_freq_axes = [
            ("NAXIS1", num_freqs as f64),
            ("NAXIS2", (num_times - 1) as f64),
            ("CRVAL1", 150e6),
            ("CDELT1", 40e3),
            ("CRPIX1", 1.0),
            ("CRVAL2", 1e9 + 1.0),
            ("CDELT2", 2.0),
            ("CRPIX2", 1.0),
        ];
        let time_freq_units = [
            ("CTYPE1", "FREQ"),
            ("CUNIT1", "Hz"),
            ("CTYPE2", "TIME"),
            ("CUNIT2", "s"),
        ];
        for (p, pol_name) in JONES_POL_NAMES.into_iter().enumerate() {
            let extname = format!("SSINS_POL={pol_name}");
            let image = read_image(&mut fptr, &extname);
            assert_eq!(image.shape(), [num_times - 1, num_freqs]);
            assert_close(
                ssins.zscore.slice(s![.., .., p]).iter(),
                image.iter().copied(),
                0.0,
            );
            assert_f64_keys(&mut fptr, &extname, &time_freq_axes);
            assert_f64_keys(&mut fptr, &extname, &[("N_BL", 3.0)]);
            assert_str_keys(&mut fptr, &extname, &time_freq_units);
            assert_str_keys(
                &mut fptr,
                &extname,
                &[
                    ("POL", pol_name),
                    ("TELESCOP", "MWA"),
                    ("INSTRUME", "SSINS"),
                    ("ORIGIN", "Birli"),
                ],
            );
        }

        let image = read_image(&mut fptr, "SSINS_FLAGS");
        assert_eq!(image.shape(), [num_times - 1, num_freqs]);
        for ((t, f), &flag) in ssins.flag_array.indexed_iter() {
            assert_eq!(image[[t, f]].to_bits(), f64::from(u8::from(flag)).to_bits());
        }
        assert!(image[[2, 4]] > 0.5 && image.sum() < 1.5);
        assert_f64_keys(&mut fptr, "SSINS_FLAGS", &time_freq_axes);
        assert_str_keys(&mut fptr, "SSINS_FLAGS", &time_freq_units);

        let image = read_image(&mut fptr, "SSINS_N_BL");
        assert_eq!(image.shape(), [num_times - 1, num_freqs]);
        for ((t, f), &num_baselines) in ssins.num_baselines_tf.indexed_iter() {
            assert_eq!(image[[t, f]].to_bits(), f64::from(num_baselines).to_bits());
        }
        // the flagged visibility is in two differences
        assert_eq!(image[[0, 2]].to_bits(), 2.0_f64.to_bits());
        assert_eq!(image[[1, 2]].to_bits(), 2.0_f64.to_bits());
        assert_eq!(image[[2, 2]].to_bits(), 3.0_f64.to_bits());
        assert_f64_keys(&mut fptr, "SSINS_N_BL", &time_freq_axes);
        assert_str_keys(&mut fptr, "SSINS_N_BL", &time_freq_units);

        // polarisation is FITS axis 1, and frequency is FITS axis 2.
        let extname = "SSINS_DIFF_MEAN_AMP_FP";
        let image = read_image(&mut fptr, extname);
        assert_eq!(image.shape(), [num_freqs, 4]);
        assert_close(ssins.diff_mean_amp_fp.iter(), image.iter().copied(), 0.0);
        assert_str_keys(&mut fptr, extname, &[("CTYPE2", "FREQ"), ("CUNIT2", "Hz")]);
        assert_f64_keys(
            &mut fptr,
            extname,
            &[
                ("NAXIS1", 4.0),
                ("NAXIS2", num_freqs as f64),
                ("CRVAL2", 150e6),
                ("CDELT2", 40e3),
                ("CRPIX2", 1.0),
            ],
        );
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod eavils_tests {
    use super::{
        test_common::{
            assert_close, assert_f64_keys, assert_str_keys, hdu_names, metrics_context, read_image,
            save_and_open, synthetic_vis,
        },
        EAVILS, JONES_POL_NAMES,
    };
    use crate::{
        marlu::{mwalib::CorrelatorContext, ndarray::s, VisSelection},
        FlagContext,
    };

    /// The values match the definition, with a flagged timestep and channel.
    #[test]
    fn test_eavils_values() {
        let (num_times, num_freqs) = (5, 4);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 2]);
        let num_baselines = ctx.antenna_pairs.len();
        let vis = synthetic_vis(&ctx);
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[1] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[2] = true;

        let eavils = EAVILS::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags);

        assert_eq!(eavils.zscore.dim(), (num_times, num_freqs, 4));
        assert_eq!(eavils.mean_amp_fp.dim(), (num_freqs, 4));
        assert_eq!(eavils.sqrt_mean_var_amp_fp.dim(), (num_freqs, 4));
        assert!(eavils.flag_array.iter().all(|&flag| !flag));
        assert_eq!(eavils.flag_array.dim(), (num_times, num_freqs));
        assert_eq!(eavils.start_time_gps_s.to_bits(), 1e9_f64.to_bits());
        assert_eq!(eavils.integration_time_s.to_bits(), 2.0_f64.to_bits());
        assert_eq!(eavils.start_freq_hz.to_bits(), 150e6_f64.to_bits());
        assert_eq!(eavils.channel_width_hz.to_bits(), 40e3_f64.to_bits());

        let unflagged_times = [0, 2, 3, 4];
        let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
        for f in 0..num_freqs {
            for p in 0..4 {
                if chan_flags[f] {
                    assert!(eavils.mean_amp_fp[[f, p]].is_nan());
                    assert!(eavils.sqrt_mean_var_amp_fp[[f, p]].is_nan());
                    assert!(eavils.zscore.slice(s![.., f, p]).iter().all(|z| z.is_nan()));
                    continue;
                }
                let amp = |t: usize, b: usize| f64::from(vis[[t, f, b]][p].norm());
                // blmean = np.mean(amp, axis=baseline)
                let blmean =
                    |t: usize| mean(&(0..num_baselines).map(|b| amp(t, b)).collect::<Vec<_>>());
                // np.mean(blmean, axis=time)
                let mean_amp = mean(&unflagged_times.map(blmean));
                // np.sqrt(np.mean(np.var(amp, axis=time, ddof=1), axis=baseline))
                let stdv = mean(
                    &(0..num_baselines)
                        .map(|b| {
                            let amps = unflagged_times.map(|t| amp(t, b));
                            let mean_amp = mean(&amps);
                            amps.map(|amp| (amp - mean_amp).powi(2)).iter().sum::<f64>()
                                / (amps.len() - 1) as f64
                        })
                        .collect::<Vec<_>>(),
                )
                .sqrt();
                assert_close([&eavils.mean_amp_fp[[f, p]]], [mean_amp], 1e-5);
                assert_close([&eavils.sqrt_mean_var_amp_fp[[f, p]]], [stdv], 1e-4);
                let zscores = (0..num_times).map(|t| {
                    if timestep_flags[t] {
                        f64::NAN
                    } else {
                        (blmean(t) - mean_amp) / stdv * (num_baselines as f64).sqrt()
                    }
                });
                // the z-score is a small difference of f32 means
                assert_close(eavils.zscore.slice(s![.., f, p]).iter(), zscores, 2e-3);
            }
        }
    }

    /// `EAVILS::new` uses the flags, times and frequencies of the selection.
    #[test]
    fn test_eavils_new_with_selection() {
        let corr_ctx = CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &[
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
            ],
        )
        .unwrap();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.timestep_range = 1..4;
        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        flag_ctx.timestep_flags[2] = true;
        flag_ctx.coarse_chan_flags[vis_sel.coarse_chan_range.start + 1] = true;
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        crate::io::read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        let eavils = EAVILS::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        let (ctx, timestep_flags, chan_flags) =
            super::selection_metadata(&corr_ctx, &vis_sel, &flag_ctx);
        assert_eq!(timestep_flags, [false, true, false]);
        assert_eq!(chan_flags, [false, false, true, true]);
        let expected =
            EAVILS::new_from_metadata(jones_array.view(), &ctx, &timestep_flags, &chan_flags);
        assert_close(
            eavils.zscore.iter(),
            expected.zscore.iter().map(|&z| f64::from(z)),
            0.0,
        );
        assert_eq!(eavils.zscore.dim(), (3, 4, 4));
        assert!(eavils
            .zscore
            .slice(s![0, ..2, ..])
            .iter()
            .all(|z| z.is_finite()));
        assert!(eavils
            .zscore
            .slice(s![1, .., ..])
            .iter()
            .all(|z| z.is_nan()));
        assert!(eavils
            .zscore
            .slice(s![.., 2.., ..])
            .iter()
            .all(|z| z.is_nan()));
        assert_eq!(
            eavils.start_time_gps_s.to_bits(),
            (corr_ctx.timesteps[1].gps_time_ms as f64 / 1000.0).to_bits()
        );
    }

    #[test]
    fn test_eavils_save_to_fits() {
        let (num_times, num_freqs) = (3, 6);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        let mut eavils = EAVILS::new_from_metadata(vis.view(), &ctx, &[false; 3], &[false; 6]);
        eavils.flag_array[[1, 5]] = true;

        let (_tmp_file, mut fptr) = save_and_open(|fptr| eavils.save_to_fits(fptr));

        assert_eq!(
            hdu_names(&mut fptr),
            [
                "EAVILS_POL=XX",
                "EAVILS_POL=XY",
                "EAVILS_POL=YX",
                "EAVILS_POL=YY",
                "EAVILS_FLAGS",
                "EAVILS_MEAN_AMP_FP",
                "EAVILS_SQRT_MEAN_VAR_AMP_FP",
            ]
        );

        let time_freq_axes = [
            ("NAXIS1", num_freqs as f64),
            ("NAXIS2", num_times as f64),
            ("CRVAL1", 150e6),
            ("CDELT1", 40e3),
            ("CRPIX1", 1.0),
            ("CRVAL2", 1e9),
            ("CDELT2", 2.0),
            ("CRPIX2", 1.0),
        ];
        let time_freq_units = [
            ("CTYPE1", "FREQ"),
            ("CUNIT1", "Hz"),
            ("CTYPE2", "TIME"),
            ("CUNIT2", "s"),
            ("TELESCOP", "MWA"),
            ("INSTRUME", "EAVILS"),
            ("ORIGIN", "Birli"),
        ];
        for (p, pol_name) in JONES_POL_NAMES.into_iter().enumerate() {
            let extname = format!("EAVILS_POL={pol_name}");
            let image = read_image(&mut fptr, &extname);
            assert_eq!(image.shape(), [num_times, num_freqs]);
            assert_close(
                eavils.zscore.slice(s![.., .., p]).iter(),
                image.iter().copied(),
                0.0,
            );
            assert_f64_keys(&mut fptr, &extname, &time_freq_axes);
            assert_str_keys(&mut fptr, &extname, &time_freq_units);
            assert_str_keys(&mut fptr, &extname, &[("POL", pol_name)]);
        }

        let image = read_image(&mut fptr, "EAVILS_FLAGS");
        assert_eq!(image.shape(), [num_times, num_freqs]);
        assert!(image[[1, 5]] > 0.5 && image.sum() < 1.5);
        assert_f64_keys(&mut fptr, "EAVILS_FLAGS", &time_freq_axes);
        assert_str_keys(&mut fptr, "EAVILS_FLAGS", &time_freq_units);

        // polarisation is FITS axis 1, and frequency is FITS axis 2.
        for (extname, array_fp) in [
            ("EAVILS_MEAN_AMP_FP", &eavils.mean_amp_fp),
            ("EAVILS_SQRT_MEAN_VAR_AMP_FP", &eavils.sqrt_mean_var_amp_fp),
        ] {
            let image = read_image(&mut fptr, extname);
            assert_eq!(image.shape(), [num_freqs, 4]);
            assert_close(array_fp.iter(), image.iter().copied(), 0.0);
            assert_str_keys(&mut fptr, extname, &[("CTYPE2", "FREQ"), ("CUNIT2", "Hz")]);
            assert_f64_keys(
                &mut fptr,
                extname,
                &[
                    ("NAXIS1", 4.0),
                    ("NAXIS2", num_freqs as f64),
                    ("CRVAL2", 150e6),
                    ("CDELT2", 40e3),
                    ("CRPIX2", 1.0),
                ],
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod flag_and_cross_tests {
    #[cfg(feature = "aoflagger")]
    use super::AOFlagMetrics;
    // some of these are only used by the tests which need the aoflagger feature
    #[cfg_attr(not(feature = "aoflagger"), allow(unused_imports))]
    use super::{
        test_common::{
            assert_close, assert_f64_keys, assert_str_keys, hdu_names, metrics_context, read_image,
            save_and_open, synthetic_vis,
        },
        CrossMetrics, DELAY_TRANSFORM_CONFIG, SSINS,
    };
    #[cfg_attr(not(feature = "aoflagger"), allow(unused_imports))]
    use crate::{
        delay_transform::delay_transform,
        marlu::{
            mwalib::CorrelatorContext,
            ndarray::{s, Array1, Array2, Array3, Axis},
            num_complex::Complex,
            Jones, VisSelection,
        },
        FlagContext,
    };
    use std::f64::consts::TAU;

    fn get_mwax_context() -> CorrelatorContext {
        CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &[
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
                "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
            ],
        )
        .unwrap()
    }

    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_crosses_only_occupancy_helper() {
        // Two timesteps, two freqs, three baselines: (0,0) auto, (0,1) cross, (1,1) auto
        let mut flags = Array3::<bool>::from_elem((2, 2, 3), false);
        // autos are flagged everywhere
        flags.slice_mut(s![.., .., 0]).fill(true);
        flags.slice_mut(s![.., .., 2]).fill(true);
        // flag cross at (t=1,f=0)
        flags[[1, 0, 1]] = true;

        let occ = AOFlagMetrics::compute_cross_occupancy(
            flags.view(),
            &[false, false],
            &[false, false],
            &[(0, 0), (0, 1), (1, 1)],
        );

        // fraction over crosses only. Only one cross baseline exists.
        assert_eq!(occ.dim(), (2, 2));
        assert_close(occ.mapv(|x| x as f32).iter(), [0.0, 0.0, 1.0, 0.0], 1e-12);
    }

    /// The occupancy is NaN if there are no cross-correlations.
    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_occupancy_without_crosses() {
        let flags = Array3::<bool>::from_elem((2, 2, 2), true);
        let occ = AOFlagMetrics::compute_cross_occupancy(
            flags.view(),
            &[false, false],
            &[false, false],
            &[(0, 0), (1, 1)],
        );
        assert!(occ.iter().all(|x| x.is_nan()));
    }

    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_aoflagmetrics_values_and_fits() {
        let (num_times, num_freqs) = (4, 5);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 2, 3]);
        let num_baselines = ctx.antenna_pairs.len();
        let cross: Vec<usize> = (0..num_baselines)
            .filter(|&b| ctx.antenna_pairs[b].0 != ctx.antenna_pairs[b].1)
            .collect();
        assert_eq!(cross.len(), 6);
        let mut flags = Array3::from_elem((num_times, num_freqs, num_baselines), false);
        // half of the crosses at (0, 0), all of them at (2, 3), and an auto
        // at (1, 1)
        for &b in &cross[..3] {
            flags[[0, 0, b]] = true;
        }
        flags.slice_mut(s![2, 3, ..]).fill(true);
        flags[[1, 1, 0]] = true;
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[3] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[4] = true;

        let metrics =
            AOFlagMetrics::new_from_metadata(flags.view(), &ctx, &timestep_flags, &chan_flags);

        for ((t, f), &occupancy) in metrics.occupancy_tf.indexed_iter() {
            let expected = match (t, f) {
                (3, _) | (_, 4) => f64::NAN,
                (0, 0) => 0.5,
                (2, 3) => 1.0,
                _ => 0.0,
            };
            assert!(
                (occupancy.is_nan() && expected.is_nan()) || (occupancy - expected).abs() < 1e-12,
                "t={t} f={f} {occupancy} != {expected}"
            );
        }

        let (_tmp_file, mut fptr) = save_and_open(|fptr| metrics.save_to_fits(fptr));
        assert_eq!(hdu_names(&mut fptr), ["AO_FLAG_METRICS"]);
        let image = read_image(&mut fptr, "AO_FLAG_METRICS");
        assert_eq!(image.shape(), [num_times, num_freqs]);
        assert_close(
            metrics.occupancy_tf.mapv(|x| x as f32).iter(),
            image.iter().copied(),
            1e-6,
        );
        assert_str_keys(
            &mut fptr,
            "AO_FLAG_METRICS",
            &[
                ("CTYPE1", "FREQ"),
                ("CUNIT1", "Hz"),
                ("CTYPE2", "TIME"),
                ("CUNIT2", "s"),
            ],
        );
        assert_f64_keys(
            &mut fptr,
            "AO_FLAG_METRICS",
            &[
                ("CRVAL1", 150e6),
                ("CDELT1", 40e3),
                ("CRPIX1", 1.0),
                ("CRVAL2", 1e9),
                ("CDELT2", 2.0),
                ("CRPIX2", 1.0),
            ],
        );
    }

    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_aoflagmetrics_new_with_selection() {
        let corr_ctx = get_mwax_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        flag_ctx.timestep_flags[0] = true;
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        // baselines are (0, 0), (0, 1), (1, 1)
        flag_array[[2, 1, 1]] = true;

        let metrics = AOFlagMetrics::new(flag_array.view(), &corr_ctx, &vis_sel, &flag_ctx);

        assert_eq!(metrics.occupancy_tf.dim(), (4, 4));
        assert!(metrics.occupancy_tf.row(0).iter().all(|x| x.is_nan()));
        assert!((metrics.occupancy_tf[[2, 1]] - 1.0).abs() < 1e-12);
        assert!(metrics.occupancy_tf[[1, 1]].abs() < 1e-12);
        assert_eq!(
            metrics.start_time_gps_s.to_bits(),
            (corr_ctx.timesteps[0].gps_time_ms as f64 / 1000.0).to_bits()
        );
        let first_chan_idx = vis_sel.coarse_chan_range.start * fine_chans_per_coarse;
        assert_eq!(
            metrics.start_freq_hz.to_bits(),
            corr_ctx.metafits_context.metafits_fine_chan_freqs_hz[first_chan_idx].to_bits()
        );
    }

    /// Visibilities whose amplitude has a ripple in frequency with a different
    /// delay for each polarisation: (none, 400 ns, 800 ns, 1200 ns).
    fn ripple_vis(num_times: usize, freqs_hz: &[f64], num_baselines: usize) -> Array3<Jones<f32>> {
        Array3::from_shape_fn((num_times, freqs_hz.len(), num_baselines), |(t, f, b)| {
            Jones::from([0, 1, 2, 3].map(|p| {
                let ripple = (TAU * freqs_hz[f] * 400e-9 * p as f64).cos();
                let amp = (1 + b) as f64 * (2.0 + 0.5 * ripple);
                // the phase changes with time, so that the coherent average
                // in time matters
                Complex::from_polar(amp as f32, 0.2 * t as f32)
            }))
        })
    }

    #[test]
    fn test_crossmetrics_values_and_fits() {
        let (num_times, num_freqs) = (4, 256);
        // tiles 10 m apart: only neighbours are within 15 m.
        let mut ctx = metrics_context(num_times, num_freqs, &[5, 1, 2, 3]);
        // tile 5 is 20 m from tile 3, so it has no short baselines.
        let num_baselines = ctx.antenna_pairs.len();
        let mut vis = ripple_vis(num_times, &ctx.fine_chan_freqs_hz, num_baselines);
        // YY of tile 1's only short baseline, (1, 2), is all zero, and
        // baseline (2, 3) has a NaN
        let bl_12 = ctx.antenna_pairs.iter().position(|&p| p == (1, 2)).unwrap();
        let bl_23 = ctx.antenna_pairs.iter().position(|&p| p == (2, 3)).unwrap();
        for t in 0..num_times {
            for f in 0..num_freqs {
                vis[[t, f, bl_12]][3] = Complex::new(0.0, 0.0);
            }
        }
        vis[[1, 7, bl_23]][1] = Complex::new(f32::NAN, 0.0);
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[2] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[9] = true;

        let cross =
            CrossMetrics::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags, 15.0);

        assert_eq!(
            cross.antenna_names,
            ["Tile001", "Tile002", "Tile003", "Tile005"]
        );
        assert_eq!(cross.antenna_ids, [1, 2, 3, 5]);
        assert_eq!(cross.antenna_nums, [1001, 1002, 1003, 1005]);
        assert_eq!(cross.start_freq_hz.to_bits(), 150e6_f64.to_bits());
        assert_eq!(cross.channel_width_hz.to_bits(), 40e3_f64.to_bits());
        assert_eq!(cross.start_time_gps_s.to_bits(), 1e9_f64.to_bits());
        assert_eq!(cross.integration_time_s.to_bits(), 2.0_f64.to_bits());

        // the delay spectrum of a baseline: the delay transform of the
        // amplitude of the mean over unflagged timesteps of its visibilities.
        let freqs_hz = Array1::from(ctx.fine_chan_freqs_hz.clone());
        let baseline_delay_spectrum = |bl_idx: usize, p: usize| -> Array1<f64> {
            let spectrum: Array1<f64> = (0..num_freqs)
                .map(|f| {
                    if chan_flags[f] {
                        return 0.0;
                    }
                    let finite: Vec<Complex<f64>> = (0..num_times)
                        .filter(|&t| !timestep_flags[t])
                        .map(|t| vis[[t, f, bl_idx]][p])
                        .filter(|v| v.re.is_finite() && v.im.is_finite())
                        .map(|v| Complex::new(f64::from(v.re), f64::from(v.im)))
                        .collect();
                    (finite.iter().sum::<Complex<f64>>() / finite.len() as f64).norm()
                })
                .collect();
            delay_transform(
                &spectrum.insert_axis(Axis(0)),
                &freqs_hz,
                &DELAY_TRANSFORM_CONFIG,
            )
            .unwrap()
            .delay_spectrum
            .row(0)
            .to_owned()
        };
        let num_delays = cross.short_ant_delay_pol_adp.dim().1;
        assert_eq!(cross.short_ant_delay_pol_adp.dim(), (4, num_delays, 4));
        for p in 0..4 {
            let delay_12 = baseline_delay_spectrum(bl_12, p);
            let delay_23 = baseline_delay_spectrum(bl_23, p);
            let expected: [Array1<f64>; 4] = if p == 3 {
                // baseline (1, 2) has no YY
                [
                    Array1::from_elem(num_delays, f64::NAN),
                    delay_23.clone(),
                    delay_23,
                    Array1::from_elem(num_delays, f64::NAN),
                ]
            } else {
                [
                    delay_12.clone(),
                    (&delay_12 + &delay_23) / 2.0,
                    delay_23,
                    Array1::from_elem(num_delays, f64::NAN),
                ]
            };
            for (a, expected) in expected.iter().enumerate() {
                assert_close(
                    cross.short_ant_delay_pol_adp.slice(s![a, .., p]).iter(),
                    expected.iter().copied(),
                    1e-5,
                );
            }
        }

        // the ripple of each polarisation peaks at its own delay
        let peak_delay_ns = |p: usize| {
            let (idx, _) = cross
                .short_ant_delay_pol_adp
                .slice(s![2, .., p])
                .iter()
                .enumerate()
                .filter(|&(i, _)| {
                    let delay_ns = cross.start_delay_ns + i as f64 * cross.delay_res_ns;
                    (delay_ns - 400.0 * p as f64).abs() < 200.0
                })
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap();
            cross.start_delay_ns + idx as f64 * cross.delay_res_ns
        };
        // XY and YX have ripples of 400 and 800 ns
        let (peak_xy_ns, peak_yx_ns) = (peak_delay_ns(1), peak_delay_ns(2));
        assert!((peak_xy_ns - 400.0).abs() < 150.0, "{peak_xy_ns}");
        assert!((peak_yx_ns - 800.0).abs() < 150.0, "{peak_yx_ns}");

        // a longer cutoff includes more baselines
        ctx.antennas[3].north_m = 40.0;
        let cross_long =
            CrossMetrics::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags, 100.0);
        assert!(cross_long
            .short_ant_delay_pol_adp
            .slice(s![.., .., ..3])
            .iter()
            .all(|x| x.is_finite()));

        let (_tmp_file, mut fptr) = save_and_open(|fptr| cross.save_to_fits(fptr));
        assert_eq!(
            hdu_names(&mut fptr),
            [
                "CROSS_DELAY_POL=XX",
                "CROSS_DELAY_POL=YY",
                "CROSS_DELAY_POL=XY",
                "CROSS_DELAY_POL=YX",
            ]
        );
        // the HDUs are in FITS order, but the array is in Jones order
        for (pol_name, p) in [("XX", 0), ("YY", 3), ("XY", 1), ("YX", 2)] {
            let extname = format!("CROSS_DELAY_POL={pol_name}");
            let image = read_image(&mut fptr, &extname);
            assert_eq!(image.shape(), [4, num_delays]);
            assert_close(
                cross.short_ant_delay_pol_adp.slice(s![.., .., p]).iter(),
                image.iter().copied(),
                0.0,
            );
            assert_str_keys(
                &mut fptr,
                &extname,
                &[
                    ("CTYPE1", "DELAY"),
                    ("CUNIT1", "ns"),
                    ("CTYPE2", "ANTENNA"),
                    ("POL", pol_name),
                    ("TELESCOP", "MWA"),
                    ("INSTRUME", "CROSS_METRICS"),
                    ("ORIGIN", "Birli"),
                ],
            );
            assert_f64_keys(
                &mut fptr,
                &extname,
                &[
                    ("CRVAL1", cross.start_delay_ns),
                    ("CRPIX1", 1.0),
                    ("DELAY_RES", cross.delay_res_ns),
                    ("N_ANTS", 4.0),
                ],
            );
        }
    }

    #[test]
    fn test_crossmetrics_new_with_selection() {
        let corr_ctx = get_mwax_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        crate::io::read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        // the two tiles are further apart than 1 m, and closer than 10 km
        let cross = CrossMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx, 1.0);
        assert_eq!(cross.short_ant_delay_pol_adp.dim().0, 2);
        assert!(cross.short_ant_delay_pol_adp.iter().all(|x| x.is_nan()));

        let cross = CrossMetrics::new(jones_array.view(), &corr_ctx, &vis_sel, &flag_ctx, 1e4);
        assert!(cross.short_ant_delay_pol_adp.iter().all(|x| x.is_finite()));
        // each tile has the same, single baseline
        assert_eq!(
            cross.short_ant_delay_pol_adp.slice(s![0, .., ..]),
            cross.short_ant_delay_pol_adp.slice(s![1, .., ..])
        );
        let meta_ants = &corr_ctx.metafits_context.antennas;
        assert_eq!(
            cross.antenna_names,
            [
                meta_ants[0].tile_name.clone(),
                meta_ants[1].tile_name.clone()
            ]
        );
        assert_eq!(cross.antenna_ids, [0, 1]);
        assert_eq!(
            cross.antenna_nums,
            [meta_ants[0].tile_id, meta_ants[1].tile_id]
        );
    }

    /// `AOFlagger` flags a strong outlier in the z-scores.
    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_ssins_flag() {
        let (num_times, num_freqs) = (33, 16);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 2]);
        let vis = synthetic_vis(&ctx);
        let flags = Array3::from_elem(vis.dim(), false);
        let mut ssins = SSINS::new_from_metadata(
            vis.view(),
            flags.view(),
            &ctx,
            &vec![false; num_times],
            &vec![false; num_freqs],
        );
        // noise-like z-scores, with one very strong outlier
        ssins.zscore = Array3::from_shape_fn((num_times - 1, num_freqs, 4), |(t, f, p)| {
            ((t * 13 + f * 7 + p * 3) % 5) as f32 * 0.1 - 0.2
        });
        ssins.zscore.slice_mut(s![20, 5, ..]).fill(1e4);
        assert!(ssins.flag_array.iter().all(|&flag| !flag));

        ssins.flag(None);

        assert_eq!(ssins.flag_array.dim(), (num_times - 1, num_freqs));
        assert!(ssins.flag_array[[20, 5]]);
        let num_flagged = ssins.flag_array.iter().filter(|&&flag| flag).count();
        assert!(num_flagged < ssins.flag_array.len() / 2, "{num_flagged}");

        // a strategy can be given by file name
        let strategy = unsafe { aoflagger_sys::cxx_aoflagger_new() }
            .FindStrategyFileGeneric(&String::from("minimal"));
        let flag_array: Array2<bool> = ssins.flag_array.clone();
        ssins.flag_array.fill(false);
        ssins.flag(Some(strategy));
        assert_eq!(ssins.flag_array, flag_array);
    }
}

/// The builders give the same metrics however the visibilities are chunked.
#[cfg(test)]
mod builder_tests {
    use super::{test_common::*, *};
    use crate::marlu::ndarray::{Array, Dimension};

    /// Split into chunks of 1, 3, then the rest of the timesteps.
    fn chunks<T>(array: &Array3<T>) -> [ArrayView3<'_, T>; 3] {
        [
            array.slice(s![0..1, .., ..]),
            array.slice(s![1..4, .., ..]),
            array.slice(s![4.., .., ..]),
        ]
    }

    /// Metadata, visibilities, flags, timestep flags and channel flags, with
    /// some of everything flagged, and one non-finite visibility.
    #[allow(clippy::type_complexity)]
    fn setup() -> (
        MetricsContext,
        Array3<Jones<f32>>,
        Array3<bool>,
        Vec<bool>,
        Vec<bool>,
    ) {
        let (num_times, num_freqs) = (7, 16);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1, 3, 5]);
        let mut vis = synthetic_vis(&ctx);
        vis[[5, 2, 1]] = Jones::nan();
        let mut flags = Array3::from_elem(vis.dim(), false);
        flags[[0, 0, 1]] = true;
        flags[[4, 3, 2]] = true;
        flags[[6, 5, 4]] = true;
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[2] = true;
        let mut chan_flags = vec![false; num_freqs];
        chan_flags[7] = true;
        (ctx, vis, flags, timestep_flags, chan_flags)
    }

    /// Bitwise equality, so that NaNs are equal.
    fn assert_same<D: Dimension>(chunked: &Array<f32, D>, whole: &Array<f32, D>) {
        assert_eq!(chunked.dim(), whole.dim());
        for (a, b) in chunked.iter().zip(whole) {
            assert_eq!(a.to_bits(), b.to_bits(), "{a} != {b}");
        }
    }

    #[test]
    fn test_autometrics_builder_chunks() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let whole = AutoMetrics::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags);

        let mut builder = AutoMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags);
        for chunk in chunks(&vis) {
            builder.push(chunk);
        }
        let chunked = builder.finish();

        assert_same(&chunked.auto_sub_aptf, &whole.auto_sub_aptf);
        assert_same(&chunked.auto_spectrum_afp, &whole.auto_spectrum_afp);
        assert_same(&chunked.auto_coeffs_apo, &whole.auto_coeffs_apo);
        assert_same(&chunked.auto_delay_afp, &whole.auto_delay_afp);
        assert_eq!(chunked.antenna_names, whole.antenna_names);
        assert_eq!(
            chunked.start_time_gps_s.to_bits(),
            whole.start_time_gps_s.to_bits()
        );
        // there is something to compare
        assert!(whole.auto_sub_aptf.iter().any(|x| x.is_finite()));
        assert!(whole.auto_sub_aptf.iter().any(|x| x.is_nan()));
    }

    #[test]
    fn test_ssins_builder_chunks() {
        let (ctx, vis, flags, timestep_flags, chan_flags) = setup();
        let whole =
            SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &timestep_flags, &chan_flags);

        let mut builder = SSINSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        for (vis_chunk, flag_chunk) in chunks(&vis).into_iter().zip(chunks(&flags)) {
            builder.push(vis_chunk, flag_chunk);
        }
        // an empty chunk changes nothing
        builder.push(vis.slice(s![0..0, .., ..]), flags.slice(s![0..0, .., ..]));
        let chunked = builder.finish();

        // the differences across the chunk boundaries are included
        assert_eq!(chunked.zscore.dim(), (6, 16, 4));
        assert_same(&chunked.zscore, &whole.zscore);
        assert_same(&chunked.diff_mean_amp_fp, &whole.diff_mean_amp_fp);
        assert_eq!(chunked.num_baselines_tf, whole.num_baselines_tf);
        assert_eq!(chunked.num_baselines, whole.num_baselines);
        assert_eq!(
            chunked.start_time_gps_s.to_bits(),
            whole.start_time_gps_s.to_bits()
        );
        // differences 0 and 3 are across chunk boundaries
        for t in [0, 3] {
            assert!(whole
                .zscore
                .slice(s![t, 0, ..])
                .iter()
                .all(|x| x.is_finite()));
            assert!(whole.num_baselines_tf.row(t).iter().any(|&n| n > 0));
        }
    }

    #[test]
    fn test_eavils_builder_chunks() {
        let (ctx, mut vis, _, timestep_flags, chan_flags) = setup();
        // EAVILS does not skip non-finite visibilities
        vis[[5, 2, 1]] = Jones::identity();
        let whole = EAVILS::new_from_metadata(vis.view(), &ctx, &timestep_flags, &chan_flags);

        let mut builder = EAVILSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        for chunk in chunks(&vis) {
            builder.push(chunk);
        }
        let chunked = builder.finish();

        assert_same(&chunked.zscore, &whole.zscore);
        assert_same(&chunked.mean_amp_fp, &whole.mean_amp_fp);
        assert_same(&chunked.sqrt_mean_var_amp_fp, &whole.sqrt_mean_var_amp_fp);
        assert!(whole.zscore.iter().any(|x| x.is_finite()));
    }

    /// The variance that is accumulated a timestep at a time is the
    /// sample variance over time of the amplitude of each baseline.
    #[test]
    fn test_eavils_builder_variance() {
        let (num_times, num_freqs) = (5, 3);
        let ctx = metrics_context(num_times, num_freqs, &[0, 1]);
        let vis = synthetic_vis(&ctx);
        let num_baselines = ctx.antenna_pairs.len();
        let mut timestep_flags = vec![false; num_times];
        timestep_flags[1] = true;

        let eavils = EAVILS::new_from_metadata(vis.view(), &ctx, &timestep_flags, &[false; 3]);

        for f in 0..num_freqs {
            for p in 0..4 {
                let mut mean_var = 0.0;
                for b in 0..num_baselines {
                    let amps: Vec<f64> = (0..num_times)
                        .filter(|&t| !timestep_flags[t])
                        .map(|t| f64::from(vis[[t, f, b]][p].norm()))
                        .collect();
                    let mean = amps.iter().sum::<f64>() / amps.len() as f64;
                    mean_var += amps.iter().map(|amp| (amp - mean).powi(2)).sum::<f64>()
                        / (amps.len() - 1) as f64
                        / num_baselines as f64;
                }
                let expected = mean_var.sqrt();
                let actual = f64::from(eavils.sqrt_mean_var_amp_fp[[f, p]]);
                assert!(
                    (actual - expected).abs() <= 1e-5 * expected,
                    "f={f} p={p}: {actual} != {expected}"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "too many timesteps")]
    fn test_autometrics_builder_too_many_timesteps() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let mut builder = AutoMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.view());
        builder.push(vis.slice(s![0..1, .., ..]));
    }

    #[test]
    #[should_panic(expected = "not all timesteps were pushed")]
    fn test_autometrics_builder_too_few_timesteps() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let mut builder = AutoMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.slice(s![0..1, .., ..]));
        builder.finish();
    }

    #[test]
    #[should_panic(expected = "too many timesteps")]
    fn test_ssins_builder_too_many_timesteps() {
        let (ctx, vis, flags, timestep_flags, chan_flags) = setup();
        let mut builder = SSINSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.view(), flags.view());
        builder.push(vis.slice(s![0..1, .., ..]), flags.slice(s![0..1, .., ..]));
    }

    #[test]
    #[should_panic(expected = "not all timesteps were pushed")]
    fn test_ssins_builder_too_few_timesteps() {
        let (ctx, vis, flags, timestep_flags, chan_flags) = setup();
        let mut builder = SSINSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.slice(s![0..1, .., ..]), flags.slice(s![0..1, .., ..]));
        builder.finish();
    }

    #[test]
    #[should_panic(expected = "too many timesteps")]
    fn test_eavils_builder_too_many_timesteps() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let mut builder = EAVILSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.view());
        builder.push(vis.slice(s![0..1, .., ..]));
    }

    #[test]
    #[should_panic(expected = "not all timesteps were pushed")]
    fn test_eavils_builder_too_few_timesteps() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let mut builder = EAVILSBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(vis.slice(s![0..1, .., ..]));
        builder.finish();
    }

    #[cfg(feature = "aoflagger")]
    #[test]
    fn test_aoflagmetrics_builder_chunks() {
        let (ctx, _, flags, timestep_flags, chan_flags) = setup();
        let whole =
            AOFlagMetrics::new_from_metadata(flags.view(), &ctx, &timestep_flags, &chan_flags);

        let mut builder = AOFlagMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags);
        for chunk in chunks(&flags) {
            builder.push(chunk);
        }
        let chunked = builder.finish();

        assert_same(
            &chunked.occupancy_tf.mapv(|x| x as f32),
            &whole.occupancy_tf.mapv(|x| x as f32),
        );
        assert!(whole.occupancy_tf.iter().any(|&x| x > 0.0));
    }

    #[cfg(feature = "aoflagger")]
    #[test]
    #[should_panic(expected = "not all timesteps were pushed")]
    fn test_aoflagmetrics_builder_too_few_timesteps() {
        let (ctx, _, flags, timestep_flags, chan_flags) = setup();
        let mut builder = AOFlagMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags);
        builder.push(flags.slice(s![0..1, .., ..]));
        builder.finish();
    }

    #[test]
    fn test_crossmetrics_builder_chunks() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        // antennas are 10 m apart, so this includes some baselines, not all
        let cutoff_m = 25.0;
        let whole = CrossMetrics::new_from_metadata(
            vis.view(),
            &ctx,
            &timestep_flags,
            &chan_flags,
            cutoff_m,
        );

        let mut builder = CrossMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags, cutoff_m);
        for chunk in chunks(&vis) {
            builder.push(chunk);
        }
        let chunked = builder.finish();

        assert_same(
            &chunked.short_ant_delay_pol_adp,
            &whole.short_ant_delay_pol_adp,
        );
        assert!(whole.short_ant_delay_pol_adp.iter().any(|x| x.is_finite()));
        assert_eq!(chunked.antenna_names, whole.antenna_names);
    }

    #[test]
    #[should_panic(expected = "not all timesteps were pushed")]
    fn test_crossmetrics_builder_too_few_timesteps() {
        let (ctx, vis, _, timestep_flags, chan_flags) = setup();
        let mut builder = CrossMetricsBuilder::new(&ctx, &timestep_flags, &chan_flags, 25.0);
        builder.push(vis.slice(s![0..1, .., ..]));
        builder.finish();
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, clippy::cognitive_complexity)]
mod ssins_tests {
    use super::{MetricsContext, SSINS};
    use crate::marlu::{
        ndarray::{s, Array3},
        num_complex::Complex,
        Jones,
    };

    /// 4 tiles with autos, in upper-triangle order.
    fn setup(num_times: usize, num_freqs: usize) -> (Array3<Jones<f32>>, MetricsContext) {
        let mut antenna_pairs = vec![];
        for i in 0..4 {
            for j in i..4 {
                antenna_pairs.push((i, j));
            }
        }
        let vis =
            Array3::from_shape_fn((num_times, num_freqs, antenna_pairs.len()), |(t, f, b)| {
                let r = ((t * 31 + f * 17 + b * 7) % 11) as f32;
                Jones::from([Complex::new(r, 1.0 - r); 4])
            });
        let ctx = MetricsContext {
            antennas: vec![],
            fine_chan_freqs_hz: (0..num_freqs).map(|f| 150e6 + 40e3 * f as f64).collect(),
            timestamps_s: (0..num_times).map(|t| 2.0 * t as f64).collect(),
            antenna_pairs,
        };
        (vis, ctx)
    }

    #[test]
    fn test_ssins_matches_reference_and_skips_autos_and_flags() {
        let (num_times, num_freqs) = (5, 3);
        let (mut vis, ctx) = setup(num_times, num_freqs);
        let mut flags = Array3::from_elem(vis.dim(), false);
        let cross: Vec<usize> = (0..ctx.antenna_pairs.len())
            .filter(|&b| ctx.antenna_pairs[b].0 != ctx.antenna_pairs[b].1)
            .collect();
        // one flagged cross visibility, and a whole flagged channel.
        flags[[2, 0, cross[0]]] = true;
        flags.slice_mut(s![.., 2, ..]).fill(true);
        let no_ts_flags = vec![false; num_times];
        let no_chan_flags = vec![false; num_freqs];
        let ssins =
            SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &no_ts_flags, &no_chan_flags);

        // reference: SSINS's INS(spectrum_type="cross", order=0).metric_ms
        let c = 4.0 / std::f64::consts::PI - 1.0;
        for f in 0..2 {
            for p in 0..4 {
                let mut metric = vec![];
                for t in 0..num_times - 1 {
                    let amps: Vec<f64> = cross
                        .iter()
                        .filter(|&&b| !flags[[t, f, b]] && !flags[[t + 1, f, b]])
                        .map(|&b| f64::from((vis[[t + 1, f, b]] - vis[[t, f, b]])[p].norm()))
                        .collect();
                    metric.push((amps.iter().sum::<f64>() / amps.len() as f64, amps.len()));
                }
                let weight: f64 = metric.iter().map(|&(_, n)| n as f64).sum();
                let mean = metric.iter().map(|&(m, n)| m * n as f64).sum::<f64>() / weight;
                for (t, &(m, n)) in metric.iter().enumerate() {
                    assert_eq!(ssins.num_baselines_tf[[t, f]] as usize, n);
                    let z = (m / mean - 1.0) * (n as f64 / c).sqrt();
                    assert!((f64::from(ssins.zscore[[t, f, p]]) - z).abs() < 1e-4);
                }
            }
        }
        // 6 cross baselines, minus the flagged one for the two differences that use it.
        assert_eq!(ssins.num_baselines, 6);
        assert_eq!(ssins.num_baselines_tf[[1, 0]], 5);
        assert_eq!(ssins.num_baselines_tf[[2, 0]], 5);
        assert_eq!(ssins.num_baselines_tf[[0, 0]], 6);
        assert!(ssins.zscore.slice(s![.., 2, ..]).iter().all(|z| z.is_nan()));
        assert!(ssins.diff_mean_amp_fp.row(2).iter().all(|z| z.is_nan()));

        // Neither the autos nor flagged visibilities change anything.
        let autos: Vec<usize> = (0..ctx.antenna_pairs.len())
            .filter(|b| !cross.contains(b))
            .collect();
        for &b in &autos {
            vis[[3, 1, b]] *= 1000.0;
        }
        vis[[2, 0, cross[0]]] *= 1000.0;
        let ssins2 =
            SSINS::new_from_metadata(vis.view(), flags.view(), &ctx, &no_ts_flags, &no_chan_flags);
        for (a, b) in ssins.zscore.iter().zip(ssins2.zscore.iter()) {
            #[allow(clippy::float_cmp)]
            let same = a == b || (a.is_nan() && b.is_nan());
            assert!(same);
        }
    }
}
