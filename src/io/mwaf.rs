//! Items related to the reading and writing of the FITS-based MWA Flag file format.
//!
//! # MWAF Format
//!
//! Similar to the GPU Fits format, mwaf files come in a set for each observation, and there is one
//! .mwaf file per gpubox (coarse channel). This file contains a binary table of all the flags for
//! that coarse channel. There is one row for each timestep-baseline combination, and there is only
//! one column. Each cell in the table contains a binary vector of flags for each fine channel in
//! the coarse channel.

use std::path::{Path, PathBuf};

use clap::crate_version;
use fitsio::{tables::ColumnDataType, tables::ColumnDescription, FitsFile};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::Itertools;
use marlu::{fitsio, fitsio_sys, mwalib, ndarray, rayon, VisSelection};
use mwalib::{
    CorrelatorContext, MWAVersion, _get_required_fits_key, _open_hdu, fits_open_hdu,
    get_required_fits_key,
};
use ndarray::prelude::*;
use rayon::prelude::*;
use regex::Regex;

use super::error::{
    IOError,
    IOError::{FitsIO, FitsOpen, InvalidFlagFilenameTemplate},
};

/// flag metadata which for a particular flag file in the set.
pub(crate) struct FlagFileHeader {
    /// The `VERSION` key from the primary hdu
    // TODO: what is this actually used for?
    pub version: String,
    /// The `OBSID` key from the primary hdu (formerly 'GPSTIME'; see
    /// <https://github.com/MWATelescope/Birli/issues/16>).
    pub obs_id: u32,
    /// The 'GPSSTART' key from the primary hdu; the GPS start time of the flags
    /// (centroid).
    pub gps_start: f64,
    /// The number of correlator fine channels per flag file, and the `NCHANS` key from the primary hdu.
    pub num_channels: u32,
    /// Total number of antennas (tiles) in the array, and the `NANTENNA` key from the primary hdu
    pub num_ants: u32,
    /// Number of timesteps in the observation, and the `NSCANS` key from the primary hdu
    pub num_timesteps: u32,
    /// The `NPOLS` key from the primary hdu
    pub num_pols: u8,
    /// The name of the software used to generate this flag file.
    pub software: String,
    /// The number of rows (timesteps × baselines), and the `NAXIS2` key from the table hdu.
    #[cfg(test)]
    pub num_rows: u32,
    /// The version of aoflagger used to generate these flags.
    pub aoflagger_version: Option<String>,
    /// The aoflagger strategy used to generate these flags.
    pub aoflagger_strategy: Option<String>,
}

/// A container for a gpubox's id, filename, channel flag count and baseline
/// flag count.
pub(crate) struct GpuboxFlags {
    /// The gpubox ID for these flags.
    pub(crate) id: usize,
    /// The mwaf filename for these flags.
    filename: PathBuf,
    channel_flag_count: Vec<u64>,
    baseline_flag_count: Vec<u64>,
}

/// A group of .mwaf files for the same observation
pub struct FlagFileSet {
    pub(crate) gpuboxes: Vec<GpuboxFlags>,
    /// Flag file header.
    pub(crate) header: FlagFileHeader,
    /// The number of baseline rows that have been written.
    row_count: u64,
    /// The number of baseline rows that are expected when writing is complete.
    expected_rows: u64,
    /// The antenna pairs of all baselines used in the flags.
    ant_pairs: Vec<(usize, usize)>,
    /// The names of the antennas used in the flags.
    ant_names: Vec<String>,
    /// The indices of the antennas used in the flags.
    ant_indices: Vec<u32>,
}

// helper to get the sorted unique antenna indices from ant pairs
pub(super) fn ant_indices(ant_pairs: &[(usize, usize)]) -> Vec<u32> {
    let unique_ant_indices = ant_pairs
        .iter()
        .flat_map(|(a, b)| [*a as u32, *b as u32].into_iter())
        .collect::<std::collections::HashSet<_>>();
    let mut unique_ant_indices_sorted = unique_ant_indices.into_iter().collect::<Vec<_>>();
    unique_ant_indices_sorted.sort_unstable();
    unique_ant_indices_sorted
}

impl FlagFileSet {
    fn get_gpubox_filenames(
        mwa_version: MWAVersion,
        filename_template: &str,
        gpubox_ids: &[usize],
    ) -> Result<Vec<GpuboxFlags>, IOError> {
        let num_percents = match mwa_version {
            MWAVersion::CorrOldLegacy | MWAVersion::CorrLegacy => 2,
            _ => 3,
        };
        let re_percents = Regex::new(format!("%{{{num_percents},}}+").as_str()).unwrap();

        if !re_percents.is_match(filename_template) {
            return Err(InvalidFlagFilenameTemplate {
                source_file: file!(),
                source_line: line!(),
                filename_template: String::from(filename_template),
            });
        }

        Ok(gpubox_ids
            .iter()
            .map(|&gpubox_id| GpuboxFlags {
                id: gpubox_id,
                filename: PathBuf::from(
                    re_percents
                        .replace(filename_template, format!("{gpubox_id:0num_percents$}"))
                        .to_string(),
                ),
                channel_flag_count: vec![],
                baseline_flag_count: vec![],
            })
            .collect())
    }

    /// Create a new set of flag files
    ///
    /// * `filename_template` is a template string which is expanded to the list
    ///   of flag files in the set, by replacing the percent (`%`) characters
    ///   with each coarse channel's zero-prefixed GPU box ID. This is to
    ///   maintain backwards compatibility with Cotter.
    ///
    ///   For MWA Ord (legacy, pre-2021) correlator observations, the GPU box ID
    ///   is the two digit correlator channel host number corresponding to
    ///   `corr_chan_number` in [`mwalib::CoarseChannel`]
    ///
    ///   For MWAX correlator observations, the GPU box ID is the three-digit
    ///   received channel number corresponding to `rec_chan_number` in
    ///   [`mwalib::CoarseChannel`].
    ///
    ///    Be sure to specify the correct number of percent characters.
    ///
    /// * `corr_ctx` is an [`mwalib::CorrelatorContext`].
    ///
    /// * `vis_sel` contains the selection of coarse channels, timesteps and
    ///    baselines corresponding to the flags that will be written.
    ///
    /// * `aoflagger_version` is the version of `aoflagger` used to generate the
    ///   flags that will be provided to this `FlagFileSet`.
    ///
    /// * `aoflagger_strategy` is the `aoflagger` strategy file used to generate
    ///   the flags that will be provided to this `FlagFileSet`.
    ///
    /// # Errors
    ///
    /// Will error with [`IOError::FitsOpen`] if there are files already present
    /// at the paths specified in filename template.
    ///
    /// Will error with [`IOError::InvalidFlagFilenameTemplate`] if an invalid
    /// flag filename template is provided (wrong number of percents).
    pub fn new(
        filename_template: &str,
        corr_ctx: &CorrelatorContext,
        vis_sel: &VisSelection,
        aoflagger_version: Option<String>,
        aoflagger_strategy: Option<String>,
    ) -> Result<Self, IOError> {
        let timestep_range = vis_sel.timestep_range.clone();
        let coarse_chan_range = vis_sel.coarse_chan_range.clone();

        let gpubox_ids = coarse_chan_range
            .into_iter()
            .map(|i| corr_ctx.coarse_chans[i].gpubox_number)
            .collect::<Vec<_>>();
        let mut gpuboxes =
            Self::get_gpubox_filenames(corr_ctx.mwa_version, filename_template, &gpubox_ids)?;

        let num_fine_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse as u32;
        let first_timestep = timestep_range.start;
        let gps_start = corr_ctx.timesteps[first_timestep].gps_time_ms as f64 / 1e3
            + corr_ctx.metafits_context.corr_int_time_ms as f64 / 2e3;

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
        let ant_indices = ant_indices(&ant_pairs);

        let num_ants = ant_indices.len();
        let num_baselines = ant_pairs.len();
        let num_rows = timestep_range.len() * num_baselines;

        let header = FlagFileHeader {
            version: "2.0".to_string(),
            obs_id: corr_ctx.metafits_context.obs_id,
            gps_start,
            num_channels: num_fine_per_coarse,
            num_ants: num_ants as u32,
            num_timesteps: timestep_range.len() as u32,
            num_pols: 1,
            // TODO: use something like https://github.com/rustyhorde/vergen
            software: format!("Birli-{}", crate_version!()),
            #[cfg(test)]
            num_rows: num_rows as u32,
            aoflagger_version,
            aoflagger_strategy,
        };

        gpuboxes.par_iter_mut().try_for_each(|gpubox| {
            // The flag counts are currently 0 capacity; make them the right
            // length.
            gpubox
                .channel_flag_count
                .resize(num_fine_per_coarse as usize, 0);
            gpubox.baseline_flag_count.resize(num_baselines, 0);

            if gpubox.filename.exists() {
                std::fs::remove_file(&gpubox.filename)?;
            }

            match FitsFile::create(&gpubox.filename).open() {
                Ok(mut fptr) => Self::write_primary_hdu(&mut fptr, &header, Some(gpubox.id)),
                Err(fits_error) => Err(FitsOpen {
                    fits_error,
                    fits_filename: gpubox.filename.clone(),
                    source_file: file!(),
                    source_line: line!(),
                }),
            }
        })?;

        let ant_names = ant_indices
            .iter()
            .map(|&i| {
                corr_ctx.metafits_context.antennas[i as usize]
                    .tile_name
                    .clone()
            })
            .collect();

        Ok(Self {
            gpuboxes,
            header,
            row_count: 0,
            expected_rows: num_rows as u64,
            ant_pairs,
            ant_names,
            ant_indices,
        })
    }

    fn write_primary_hdu(
        fptr: &mut FitsFile,
        header: &FlagFileHeader,
        gpubox_id: Option<usize>,
    ) -> Result<(), IOError> {
        let FlagFileHeader {
            version,
            obs_id,
            gps_start,
            num_channels,
            num_ants,
            num_timesteps,
            num_pols,
            software,
            #[cfg(test)]
                num_rows: _,
            aoflagger_version,
            aoflagger_strategy,
        } = header;
        let hdu = fits_open_hdu!(fptr, 0)?;

        // Signal that we're using long strings.
        let mut status = 0;
        unsafe {
            // ffplsw = fits_write_key_longwarn
            fitsio_sys::ffplsw(
                fptr.as_raw(), /* I - FITS file pointer  */
                &mut status,   /* IO - error status       */
            );
            fitsio::errors::check_status(status)?;
        }

        hdu.write_key(fptr, "VERSION", version.as_str())?;
        hdu.write_key(fptr, "OBSID", *obs_id)?;
        // For some silly reason, writing `gps_start` as a float truncates;
        // writing a string works fine.
        hdu.write_key(fptr, "GPSSTART", gps_start.to_string())?;
        hdu.write_key(fptr, "NCHANS", *num_channels)?;
        hdu.write_key(fptr, "NANTENNA", *num_ants)?;
        hdu.write_key(fptr, "NSCANS", *num_timesteps)?;
        hdu.write_key(fptr, "NPOLS", *num_pols as u32)?;
        if let Some(gpubox_id) = gpubox_id {
            hdu.write_key(fptr, "GPUBOXNO", gpubox_id as u32)?;
        }
        hdu.write_key(fptr, "SOFTWARE", software.as_str())?;
        if let Some(aoflagger_version) = aoflagger_version {
            hdu.write_key(fptr, "AO_VER", aoflagger_version.as_str())?;
        }
        if let Some(aoflagger_strategy) = aoflagger_strategy.as_deref() {
            // Write AO_STRAT as a long string.
            let key_name = std::ffi::CString::new("AO_STRAT").unwrap();
            let value = std::ffi::CString::new(aoflagger_strategy).unwrap();
            let comment = std::ffi::CString::new("Strategy file used by AOFlagger").unwrap();
            unsafe {
                // ffpkls = fits_write_key_longstr
                fitsio_sys::ffpkls(
                    fptr.as_raw(),     /* I - FITS file pointer        */
                    key_name.as_ptr(), /* I - name of keyword to write */
                    value.as_ptr(),    /* I - keyword value            */
                    comment.as_ptr(),  /* I - keyword comment          */
                    &mut status,       /* IO - error status            */
                );
                fitsio::errors::check_status(status)?;
            }
        }

        {
            // Write CMDLINE as a long string.
            let cmd_line = std::env::args().collect_vec().join(" ");
            let key_name = std::ffi::CString::new("CMDLINE").unwrap();
            let value = std::ffi::CString::new(cmd_line).unwrap();
            let comment = std::ffi::CString::new("Command-line call").unwrap();
            unsafe {
                // ffpkls = fits_write_key_longstr
                fitsio_sys::ffpkls(
                    fptr.as_raw(),     /* I - FITS file pointer        */
                    key_name.as_ptr(), /* I - name of keyword to write */
                    value.as_ptr(),    /* I - keyword value            */
                    comment.as_ptr(),  /* I - keyword comment          */
                    &mut status,       /* IO - error status            */
                );
                fitsio::errors::check_status(status)?;
            }
        }

        Ok(())
    }

    /// Write the provided flags to disk, given an ndarray of boolean flags for
    /// the observation into a file for each `gpubox_id`. Progress bars can also
    /// be drawn.
    ///
    /// The filename template should contain two or 3 percentage (`%`)
    /// characters which will be replaced by the gpubox id or channel number
    /// (depending on correlator type). See [`FlagFileSet::new`]
    ///
    /// # Errors
    ///
    /// Will error if the gpubox ids this flagset was initialized with is not
    /// contained in the provided [`mwalib::CorrelatorContext`].
    ///
    pub fn write_flag_array(
        &mut self,
        flag_array: ArrayView3<bool>,
        draw_progress: bool,
    ) -> Result<(), IOError> {
        let flag_dims = flag_array.dim();
        let num_timesteps = flag_dims.0;
        let num_baselines = flag_dims.2;
        let num_fine_chans_per_coarse = self.header.num_channels as usize;
        assert_eq!(num_fine_chans_per_coarse * self.gpuboxes.len(), flag_dims.1);

        let multi_progress = MultiProgress::with_draw_target(if draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        });
        // progress bars detailing the write progress of each mwaf file.
        let progress_bars = self
            .gpuboxes
            .iter()
            .map(|gpubox| {
                multi_progress.add({
                    ProgressBar::new((num_timesteps * num_baselines) as _)
                        .with_style(
                            ProgressStyle::default_bar()
                                .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                                .unwrap()
                                .progress_chars("=> "),
                        )
                        .with_position(0)
                        .with_message(format!("mwaf {:03}", gpubox.id))
                })
            })
            .collect::<Vec<_>>();

        self.gpuboxes
            .par_iter_mut()
            .zip(
                flag_array
                    .axis_chunks_iter(Axis(1), num_fine_chans_per_coarse)
                    .into_par_iter(),
            )
            .zip(progress_bars)
            .try_for_each(|((gpubox, flag_coarse_chan_view), channel_progress)| {
                Self::write_flag_array_inner(
                    &gpubox.filename,
                    flag_coarse_chan_view,
                    &channel_progress,
                    num_fine_chans_per_coarse,
                    &mut gpubox.channel_flag_count,
                    &mut gpubox.baseline_flag_count,
                )
            })?;

        self.row_count += (flag_array.len_of(Axis(0)) * flag_array.len_of(Axis(2))) as u64;

        Ok(())
    }

    /// This fallible function is run in parallel from `write_flag_array`.
    fn write_flag_array_inner(
        filename: &Path,
        flag_coarse_chan_view: ArrayView3<bool>,
        channel_progress: &ProgressBar,
        num_fine_chans_per_coarse: usize,
        channel_flag_counts: &mut [u64],
        baseline_flag_counts: &mut [u64],
    ) -> Result<(), IOError> {
        let mut fptr = FitsFile::edit(filename)?;
        // If the FLAGS HDU doesn't exist, create it.
        let mut row_idx = if fptr.hdu("FLAGS").is_err() {
            let col = ColumnDescription::new("FLAGS")
                .with_type(ColumnDataType::Bit)
                .that_repeats(num_fine_chans_per_coarse)
                .create()?;
            fptr.create_table("FLAGS", &[col])?;
            // Start the row index from 0.
            0
        } else {
            // Start the row index from where we left off (NAXIS2).
            let hdu1 = fits_open_hdu!(&mut fptr, 1)?;
            get_required_fits_key!(&mut fptr, &hdu1, "NAXIS2")?
        };

        let mut status = 0;

        let mut flag_cell = vec![0; flag_coarse_chan_view.len_of(Axis(1))];
        for flag_timestep_view in flag_coarse_chan_view.outer_iter() {
            for (flag_baseline_view, baseline_flag_count) in flag_timestep_view
                .axis_iter(Axis(1))
                .zip_eq(baseline_flag_counts.iter_mut())
            {
                flag_cell
                    .iter_mut()
                    .zip_eq(flag_baseline_view)
                    .zip_eq(channel_flag_counts.iter_mut())
                    .for_each(|((a, &b), count)| {
                        *a = i8::from(b);
                        if b {
                            *count += 1;
                        }
                    });
                *baseline_flag_count += flag_cell.iter().map(|i| *i as u64).sum::<u64>();

                unsafe {
                    fitsio_sys::ffpclx(
                        fptr.as_raw(),
                        1,
                        1 + row_idx as i64,
                        1,
                        flag_cell.len() as i64,
                        flag_cell.as_mut_ptr(),
                        &mut status,
                    );
                }

                fitsio::errors::check_status(status).map_err(|e| FitsIO {
                    fits_error: e,
                    fits_filename: fptr.filename.clone(),
                    hdu_num: 2,
                    source_file: file!(),
                    source_line: line!(),
                })?;

                row_idx += 1;
                channel_progress.inc(1);
            }
        }
        channel_progress.finish();

        Ok(())
    }

    /// Properly close the flag files.
    ///
    /// # Errors
    ///
    /// Will error if fewer than expected number of flag rows are written, or
    /// there are problems in creating fits tables in the mwaf files.
    ///
    pub fn finalise(self) -> Result<(), IOError> {
        // Complain if not enough rows were written.
        if self.row_count != self.expected_rows {
            return Err(IOError::MwafIncorrectFlagCount {
                count: self.row_count,
                expected: self.expected_rows,
            });
        }

        // Write occupancy stats and tile info into new HDUs.
        self.gpuboxes.into_par_iter().try_for_each(|gpubox| {
            Self::finalise_inner(
                gpubox,
                self.expected_rows,
                &self.ant_pairs,
                &self.ant_names,
                &self.ant_indices,
            )
        })?;

        Ok(())
    }

    /// This fallible function is run in parallel from `finalise`.
    fn finalise_inner(
        gpubox: GpuboxFlags,
        total_row_count: u64,
        ant_pairs: &[(usize, usize)],
        ant_names: &[String],
        ant_indices: &[u32],
    ) -> Result<(), IOError> {
        let mut fptr = FitsFile::edit(gpubox.filename)?;

        {
            let index_col = ColumnDescription::new("Index")
                .with_type(ColumnDataType::Int)
                .create()?;
            let count_col = ColumnDescription::new("Count")
                .with_type(ColumnDataType::Long)
                .create()?;
            let occ_col = ColumnDescription::new("Occupancy")
                .with_type(ColumnDataType::Double)
                .create()?;
            let hdu = fptr.create_table("CH_OCC", &[index_col, count_col, occ_col])?;
            hdu.write_col(
                &mut fptr,
                "Index",
                &(0..gpubox.channel_flag_count.len())
                    .map(|i| i as u32)
                    .collect::<Vec<u32>>(),
            )?;
            hdu.write_col(&mut fptr, "Count", &gpubox.channel_flag_count)?;
            hdu.write_col(
                &mut fptr,
                "Occupancy",
                &gpubox
                    .channel_flag_count
                    .iter()
                    .map(|&c| c as f64 / total_row_count as f64)
                    .collect::<Vec<f64>>(),
            )?;
        }

        {
            let index_col = ColumnDescription::new("Index")
                .with_type(ColumnDataType::Int)
                .create()?;
            let ant1_col = ColumnDescription::new("Antenna1")
                .with_type(ColumnDataType::Int)
                .create()?;
            let ant2_col = ColumnDescription::new("Antenna2")
                .with_type(ColumnDataType::Int)
                .create()?;
            let count_col = ColumnDescription::new("Count")
                .with_type(ColumnDataType::Long)
                .create()?;
            let occ_col = ColumnDescription::new("Occupancy")
                .with_type(ColumnDataType::Double)
                .create()?;
            let hdu = fptr.create_table(
                "BL_OCC",
                &[index_col, ant1_col, ant2_col, count_col, occ_col],
            )?;
            hdu.write_col(
                &mut fptr,
                "Index",
                &(0..ant_pairs.len()).map(|i| i as u32).collect::<Vec<u32>>(),
            )?;
            let (ant1s, ant2s): (Vec<u32>, Vec<u32>) = ant_pairs
                .iter()
                .map(|&(a1, a2)| (a1 as u32, a2 as u32))
                .unzip();
            hdu.write_col(&mut fptr, "Antenna1", &ant1s)?;
            hdu.write_col(&mut fptr, "Antenna2", &ant2s)?;
            hdu.write_col(&mut fptr, "Count", &gpubox.baseline_flag_count)?;
            let num_timesteps = total_row_count as usize / gpubox.baseline_flag_count.len();
            let num_channels = gpubox.channel_flag_count.len();
            hdu.write_col(
                &mut fptr,
                "Occupancy",
                &gpubox
                    .baseline_flag_count
                    .into_iter()
                    .map(|c| c as f64 / (num_timesteps * num_channels) as f64)
                    .collect::<Vec<f64>>(),
            )?;
        }

        {
            let index_col = ColumnDescription::new("Antenna")
                .with_type(ColumnDataType::Int)
                .create()?;
            let name_col = ColumnDescription::new("TileName")
                .with_type(ColumnDataType::String)
                .that_repeats(8)
                .create()?;
            let hdu = fptr.create_table("TILES", &[index_col, name_col])?;
            hdu.write_col(&mut fptr, "Antenna", ant_indices)?;
            hdu.write_col(&mut fptr, "TileName", ant_names)?;
        }

        Ok(())
    }

    #[cfg(test)]
    fn read_header(fptr: &mut FitsFile) -> Result<(FlagFileHeader, Option<u32>), ReadMwafError> {
        use mwalib::{_get_optional_fits_key, get_optional_fits_key};

        let hdu0 = fits_open_hdu!(fptr, 0)?;
        let version = get_required_fits_key!(fptr, &hdu0, "VERSION")?;
        let obs_id = get_required_fits_key!(fptr, &hdu0, "OBSID")?;
        let gps_start = get_required_fits_key!(fptr, &hdu0, "GPSSTART")?;
        let num_channels = get_required_fits_key!(fptr, &hdu0, "NCHANS")?;
        let num_ants = get_required_fits_key!(fptr, &hdu0, "NANTENNA")?;
        let num_timesteps = get_required_fits_key!(fptr, &hdu0, "NSCANS")?;
        let num_pols = get_required_fits_key!(fptr, &hdu0, "NPOLS")?;
        let gpubox_id = get_optional_fits_key!(fptr, &hdu0, "GPUBOXNO")?;
        let software = get_required_fits_key!(fptr, &hdu0, "SOFTWARE")?;
        let aoflagger_version = get_optional_fits_key!(fptr, &hdu0, "AO_VER")?;
        let aoflagger_strategy = get_optional_fits_key!(fptr, &hdu0, "AO_STRAT")?;

        let hdu1 = fits_open_hdu!(fptr, 1)?;
        let num_rows = get_required_fits_key!(fptr, &hdu1, "NAXIS2")?;

        let header = FlagFileHeader {
            version,
            obs_id,
            gps_start,
            num_channels,
            num_ants,
            num_timesteps,
            num_pols,
            software,
            num_rows,
            aoflagger_version,
            aoflagger_strategy,
        };
        let baselines = (header.num_ants * (header.num_ants + 1)) / 2;
        if header.num_rows != header.num_timesteps * baselines {
            return Err(ReadMwafError::Generic(format!(
                "File {:?}: Expected NSCANS * NANTENNA * (NANTENNA+1) / 2 = NAXIS2, found {} * {} != {}",
                fptr.filename,
                header.num_timesteps, baselines, header.num_rows
            )));
        }
        Ok((header, gpubox_id))
    }

    /// Open an existing set of flag files, given an observation's MWA Version,
    /// the flag filename template, and a list of gpubox ids. Only used in
    /// testing.
    #[cfg(test)]
    pub(crate) fn open(
        filename_template: &str,
        gpubox_ids: &[usize],
        mwa_version: MWAVersion,
    ) -> Result<Self, String> {
        let gpuboxes = Self::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)
            .map_err(|e| e.to_string())?;
        let mut header = None;
        for gpubox in &gpuboxes {
            match FitsFile::open(&gpubox.filename) {
                Ok(mut fptr) => {
                    if header.is_none() {
                        header = Some(Self::read_header(&mut fptr).map_err(|e| e.to_string())?.0);
                    }
                }
                Err(fits_error) => {
                    return Err(FitsOpen {
                        fits_error,
                        fits_filename: gpubox.filename.clone(),
                        source_file: file!(),
                        source_line: line!(),
                    }
                    .to_string())
                }
            }
        }

        let header = header.unwrap();
        let num_ants = header.num_ants as usize;

        Ok(Self {
            gpuboxes,
            header,
            row_count: 0,
            expected_rows: 0,
            ant_pairs: vec![],
            ant_names: vec![String::new(); num_ants],
            ant_indices: vec![0; num_ants],
        })
    }

    #[cfg(test)]
    #[allow(dead_code)]
    pub(crate) fn open_cotter(
        filename_template: &str,
        gpubox_ids: &[usize],
        mwa_version: MWAVersion,
    ) -> Result<(Self, String), String> {
        let gpuboxes = Self::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)
            .map_err(|e| e.to_string())?;
        let mut header = None;
        let mut gpubox_ids = Vec::with_capacity(gpuboxes.len());
        let mut date = None;

        for gpubox in &gpuboxes {
            match FitsFile::open(&gpubox.filename) {
                Ok(mut fptr) => {
                    let hdu0 = fits_open_hdu!(&mut fptr, 0).unwrap();
                    let version = get_required_fits_key!(&mut fptr, &hdu0, "VERSION").unwrap();
                    let obs_id = get_required_fits_key!(&mut fptr, &hdu0, "GPSTIME").unwrap();
                    let num_channels = get_required_fits_key!(&mut fptr, &hdu0, "NCHANS").unwrap();
                    let num_ants = get_required_fits_key!(&mut fptr, &hdu0, "NANTENNA").unwrap();
                    let num_timesteps = get_required_fits_key!(&mut fptr, &hdu0, "NSCANS").unwrap();
                    let num_pols = get_required_fits_key!(&mut fptr, &hdu0, "NPOLS").unwrap();
                    let gpubox_id: u32 =
                        get_required_fits_key!(&mut fptr, &hdu0, "GPUBOXNO").unwrap();
                    let software = get_required_fits_key!(&mut fptr, &hdu0, "COTVER").unwrap();
                    let fdate = get_required_fits_key!(&mut fptr, &hdu0, "COTVDATE").unwrap();

                    let hdu1 = fits_open_hdu!(&mut fptr, 1).unwrap();
                    let num_rows = get_required_fits_key!(&mut fptr, &hdu1, "NAXIS2").unwrap();

                    assert_eq!(gpubox.id, gpubox_id as usize);
                    gpubox_ids.push(gpubox_id);

                    if header.is_none() {
                        header = Some(FlagFileHeader {
                            version,
                            obs_id,
                            gps_start: obs_id as f64,
                            num_channels,
                            num_ants,
                            num_timesteps,
                            num_pols,
                            software,
                            num_rows,
                            aoflagger_version: None,
                            aoflagger_strategy: None,
                        });
                    }

                    if date.is_none() {
                        date = Some(fdate);
                    }
                }
                Err(fits_error) => {
                    return Err(FitsOpen {
                        fits_error,
                        fits_filename: gpubox.filename.clone(),
                        source_file: file!(),
                        source_line: line!(),
                    }
                    .to_string())
                }
            }
        }

        Ok((
            Self {
                gpuboxes,
                header: header.unwrap(),
                row_count: 0,
                expected_rows: 0,
                ant_pairs: vec![],
                ant_names: vec![],
                ant_indices: vec![],
            },
            date.unwrap(),
        ))
    }

    #[cfg(test)]
    pub(crate) fn read_flags(&self) -> Result<Array3<i8>, IOError> {
        let gpubox = &self.gpuboxes[0];
        let mut fptr = FitsFile::open(&gpubox.filename)?;
        let hdu = fits_open_hdu!(&mut fptr, 0)?;
        let num_timesteps = get_required_fits_key!(&mut fptr, &hdu, "NSCANS")?;
        let num_channels_per_mwaf: usize = get_required_fits_key!(&mut fptr, &hdu, "NCHANS")?;
        let total_num_channels = num_channels_per_mwaf * self.gpuboxes.len();
        let hdu = fits_open_hdu!(&mut fptr, 1)?;
        let num_rows: usize = get_required_fits_key!(&mut fptr, &hdu, "NAXIS2")?;
        assert_eq!(
            num_rows % num_timesteps,
            0,
            "num_rows={num_rows} should be a multiple of num_timesteps={num_timesteps}"
        );
        let num_baselines = num_rows / num_timesteps;
        let hdu = fits_open_hdu!(&mut fptr, 1)?;

        let mut out = Array3::zeros((num_timesteps, num_baselines, total_num_channels));
        drop(fptr);
        drop(hdu);

        let mut row_flags = Array1::zeros(num_channels_per_mwaf);
        for (i_gpubox, gpubox) in self.gpuboxes.iter().enumerate() {
            let mut fptr = FitsFile::open(&gpubox.filename)?;
            fits_open_hdu!(&mut fptr, 1)?;
            let mut status = 0;
            // cfitsio won't allow you to read everything in at once. So we read
            // row-by-row. Sad face.
            for i_row in 0..num_rows {
                unsafe {
                    fitsio_sys::ffgcx(
                        fptr.as_raw(),
                        1,
                        1 + i_row as i64,
                        1,
                        num_channels_per_mwaf as i64,
                        row_flags.as_mut_ptr(),
                        &mut status,
                    );
                }
                fitsio::errors::check_status(status).map_err(|e| FitsIO {
                    fits_error: e,
                    fits_filename: fptr.filename.clone(),
                    hdu_num: 1,
                    source_file: file!(),
                    source_line: line!(),
                })?;

                out.slice_mut(s![
                    i_row / num_baselines,
                    i_row % num_baselines,
                    i_gpubox * num_channels_per_mwaf..(i_gpubox + 1) * num_channels_per_mwaf,
                ])
                .assign(&row_flags);
            }
        }

        Ok(out)
    }

    #[cfg(test)]
    pub(crate) fn read_ch_occ(&self) -> Result<(Vec<u32>, Vec<f32>), IOError> {
        use itertools::izip;

        let gpubox = &self.gpuboxes[0];
        let mut fptr = FitsFile::open(&gpubox.filename)?;
        let hdu = fits_open_hdu!(&mut fptr, 2)?;
        let num_rows: usize = get_required_fits_key!(&mut fptr, &hdu, "NAXIS2")?;
        let total_num_channels = num_rows * self.gpuboxes.len();

        let mut out_count = vec![0_u32; total_num_channels];
        let mut out_occ = vec![0.; total_num_channels];
        drop(fptr);
        drop(hdu);

        for (gpubox, out_count, out_occ) in izip!(
            self.gpuboxes.iter(),
            out_count.chunks_mut(num_rows),
            out_occ.chunks_mut(num_rows),
        ) {
            let mut fptr = FitsFile::open(&gpubox.filename)?;
            let hdu = fits_open_hdu!(&mut fptr, 2)?;
            let tmp_count: Vec<u32> = hdu.read_col(&mut fptr, "Count").unwrap();
            assert_eq!(tmp_count.len(), out_count.len());
            out_count.copy_from_slice(&tmp_count);
            let tmp_occ: Vec<f32> = hdu.read_col(&mut fptr, "Occupancy").unwrap();
            assert_eq!(tmp_occ.len(), out_occ.len());
            out_occ.copy_from_slice(&tmp_occ);
        }

        Ok((out_count, out_occ))
    }

    #[cfg(test)]
    #[allow(clippy::type_complexity)]
    pub(crate) fn read_bl_occ(
        &self,
    ) -> Result<(Vec<(u32, u32)>, Array2<u32>, Array2<f32>), IOError> {
        use itertools::izip;

        let gpubox = &self.gpuboxes[0];
        let mut fptr = FitsFile::open(&gpubox.filename)?;
        let hdu = fits_open_hdu!(&mut fptr, 3)?;
        let num_rows: usize = get_required_fits_key!(&mut fptr, &hdu, "NAXIS2")?;

        let ant1s: Vec<u32> = hdu.read_col(&mut fptr, "Antenna1").unwrap();
        let ant2s: Vec<u32> = hdu.read_col(&mut fptr, "Antenna2").unwrap();
        let ant_pairs = izip!(ant1s.into_iter(), ant2s.into_iter()).collect::<Vec<_>>();
        let num_coarse_chans = self.gpuboxes.len();
        let total_num_baselines = num_rows * num_coarse_chans;

        let mut out_count = vec![0_u32; total_num_baselines];
        let mut out_occ = vec![0.; total_num_baselines];
        drop(fptr);
        drop(hdu);

        for (gpubox, out_count, out_occ) in izip!(
            self.gpuboxes.iter(),
            out_count.chunks_mut(num_rows),
            out_occ.chunks_mut(num_rows),
        ) {
            let mut fptr = FitsFile::open(&gpubox.filename)?;
            let hdu = fits_open_hdu!(&mut fptr, 3)?;
            let tmp_count: Vec<u32> = hdu.read_col(&mut fptr, "Count").unwrap();
            assert_eq!(tmp_count.len(), out_count.len());
            out_count.copy_from_slice(&tmp_count);
            let tmp_occ: Vec<f32> = hdu.read_col(&mut fptr, "Occupancy").unwrap();
            assert_eq!(tmp_occ.len(), out_occ.len());
            out_occ.copy_from_slice(&tmp_occ);
        }

        Ok((
            ant_pairs,
            Array2::from_shape_vec((num_coarse_chans, num_rows), out_count).unwrap(),
            Array2::from_shape_vec((num_coarse_chans, num_rows), out_occ).unwrap(),
        ))
    }
}

#[cfg(test)]
#[derive(thiserror::Error, Debug)]
pub(crate) enum ReadMwafError {
    #[error(transparent)]
    Fits(#[from] mwalib::FitsError),

    #[error("{0}")]
    Generic(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        io::error::IOError::InvalidFlagFilenameTemplate,
        test_common::{get_mwa_ord_context, get_mwax_context},
        FlagContext,
    };
    use approx::{abs_diff_eq, assert_abs_diff_eq};
    use fitsio::FitsFile;
    use itertools::izip;
    use marlu::{
        fitsio,
        mwalib::{
            _get_optional_fits_key, _get_required_fits_key, fits_open_hdu, get_optional_fits_key,
            get_required_fits_key,
        },
    };
    use std::fs::File;
    use std::path::Path;
    use tempfile::tempdir;

    #[test]
    fn test_flagfileset_enforces_percents_in_filename_template() {
        let mwax_context = get_mwax_context();
        let mwa_ord_context = get_mwa_ord_context();

        macro_rules! test_percent_enforcement {
            ($template_suffix:expr, $corr_ctx:expr, $timestep_range:expr, $coarse_chan_range:expr, $expected:pat) => {
                let vis_sel = VisSelection::from_mwalib(&$corr_ctx).unwrap();
                let tmp_dir = tempdir().unwrap();
                assert!(matches!(
                    FlagFileSet::new(
                        tmp_dir.path().join($template_suffix).to_str().unwrap(),
                        $corr_ctx,
                        &vis_sel,
                        None,
                        None
                    ),
                    $expected
                ))
            };
        }
        test_percent_enforcement!(
            "mwax_no_percents.mwaf",
            &mwax_context,
            0..1,
            0..3,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_no_percents.mwaf",
            &mwa_ord_context,
            0..1,
            1..3,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwax_insufficient_percents_2_%%.mwaf",
            &mwax_context,
            0..1,
            0..3,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_sufficient_percents_2_%%.mwaf",
            &mwa_ord_context,
            0..1,
            1..3,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            "mwax_sufficient_percents_3_%%%.mwaf",
            &mwax_context,
            0..1,
            0..3,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_sufficient_percents_3_%%%.mwaf",
            &mwax_context,
            0..1,
            0..3,
            Ok(FlagFileSet { .. })
        );
    }

    #[test]
    fn test_flagfileset_new_passes_with_existing() {
        let context = get_mwax_context();
        let mut vis_sel = VisSelection::from_mwalib(&context).unwrap();
        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        let ok_gpuboxes = 0..2;
        let colliding_gpuboxes = 1..3;

        for gpubox_id in colliding_gpuboxes.clone() {
            let colliding_filename = tmp_dir.path().join(format!("Flagfile{gpubox_id:03}.mwaf"));
            File::create(colliding_filename.to_str().unwrap()).unwrap();
        }

        vis_sel.coarse_chan_range = ok_gpuboxes;
        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &context,
                &vis_sel,
                None,
                None
            ),
            Ok(FlagFileSet { .. })
        ));
        vis_sel.coarse_chan_range = colliding_gpuboxes;
        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &context,
                &vis_sel,
                None,
                None
            ),
            Ok(FlagFileSet { .. })
        ));
    }

    #[test]
    fn test_flagfileset_open_fails_with_non_existing() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("FlagfileNONEXISTING%%%.mwaf");

        let gpuboxes = gpubox_ids[..1].to_vec();

        let result = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpuboxes,
            context.mwa_version,
        );
        assert!(result.is_err());
        assert!(result.err().unwrap().contains("Couldn't open"));
    }

    #[test]
    fn test_read_headers() {
        let test_dir = Path::new("tests/data/1247842824_flags/");

        let context = CorrelatorContext::new(
            test_dir.join("1247842824.metafits"),
            &[test_dir.join("1247842824_20190722150008_gpubox01_00.fits")],
        )
        .unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let filename_template = &test_dir.join("FlagfileCotterMWA%%.mwaf");
        let gpuboxes = FlagFileSet::get_gpubox_filenames(
            MWAVersion::CorrLegacy,
            filename_template.as_os_str().to_str().unwrap(),
            &gpubox_ids,
        )
        .unwrap();

        for gpubox in gpuboxes {
            let mut fptr = FitsFile::open(gpubox.filename).unwrap();
            let hdu0 = fptr.primary_hdu().unwrap();

            let version: String = get_required_fits_key!(&mut fptr, &hdu0, "VERSION").unwrap();
            let obs_id: Option<i32> = get_optional_fits_key!(&mut fptr, &hdu0, "OBSID").unwrap();
            let gps_start: Option<f64> =
                get_optional_fits_key!(&mut fptr, &hdu0, "GPSSTART").unwrap();
            let num_channels: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NCHANS").unwrap();
            let num_ants: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NANTENNA").unwrap();
            let num_timesteps: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NSCANS").unwrap();
            let num_pols: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NPOLS").unwrap();
            let software: Option<String> =
                get_optional_fits_key!(&mut fptr, &hdu0, "SOFTWARE").unwrap();
            let aoflagger_version: Option<String> =
                get_optional_fits_key!(&mut fptr, &hdu0, "AO_VER").unwrap();
            let aoflagger_strategy: Option<String> =
                get_optional_fits_key!(&mut fptr, &hdu0, "AO_STRAT").unwrap();

            let hdu1 = fits_open_hdu!(&mut fptr, 1).unwrap();
            let bytes_per_row: i32 = get_required_fits_key!(&mut fptr, &hdu1, "NAXIS1").unwrap();
            let num_rows: i32 = get_required_fits_key!(&mut fptr, &hdu1, "NAXIS2").unwrap();

            assert_eq!(version, "1.0");
            assert_eq!(obs_id, None);
            assert_eq!(gps_start, None);
            assert_eq!(num_channels, 128);
            assert_eq!(num_ants, 128);
            assert_eq!(num_timesteps, 2);
            assert_eq!(num_pols, 1);
            assert_eq!(software, None);
            assert_eq!(aoflagger_version, None);
            assert_eq!(aoflagger_strategy, None);
            assert_eq!(bytes_per_row, 16);
            assert_eq!(num_rows, 16512);
            // Test the other keys cotter writes but we no longer write.
            let gpstime: f64 = fptr
                .primary_hdu()
                .unwrap()
                .read_key(&mut fptr, "GPSTIME")
                .unwrap();
            assert!((gpstime - 1247842824.0).abs() < f64::EPSILON);

            let cotver: String = fptr
                .primary_hdu()
                .unwrap()
                .read_key(&mut fptr, "COTVER")
                .unwrap();
            assert_eq!(cotver, "4.5");
            let cotdate: String = fptr
                .primary_hdu()
                .unwrap()
                .read_key(&mut fptr, "COTVDATE")
                .unwrap();
            assert_eq!(cotdate, "2020-08-10");
        }
    }

    #[test]
    fn test_write_primary_hdu() {
        use std::collections::BTreeMap;

        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();
        let first_timestep = context.common_timestep_indices[0];
        let gps_start = context.timesteps[first_timestep].gps_time_ms as f64 / 1e3
            + context.metafits_context.corr_int_time_ms as f64 / 2e3;

        let tmp_dir = tempdir().unwrap();
        let mut gpubox_paths = BTreeMap::new();
        for &gpubox_id in &gpubox_ids {
            gpubox_paths.insert(
                gpubox_id,
                tmp_dir.path().join(format!("Flagfile{gpubox_id:03}.mwaf")),
            );
        }

        #[cfg(feature = "aoflagger")]
        let (aoflagger_version, aoflagger_strategy) = {
            let mut major = 0;
            let mut minor = 0;
            let mut subminor = 0;
            unsafe {
                aoflagger_sys::cxx_aoflagger_new().GetVersion(
                    &mut major,
                    &mut minor,
                    &mut subminor,
                );
            }
            (
                Some(format!("v{major}.{minor}.{subminor}")),
                Some("amazing_flagging_strategy.lua".to_string()),
            )
        };
        #[cfg(not(feature = "aoflagger"))]
        let (aoflagger_version, aoflagger_strategy) = (None, None);

        let num_fine_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse as u32;

        {
            for (&gpubox_id, path) in &gpubox_paths {
                let mut fptr = FitsFile::create(path).open().unwrap();
                let header = FlagFileHeader {
                    version: "2.0".to_string(),
                    obs_id: context.metafits_context.obs_id,
                    gps_start,
                    num_channels: num_fine_per_coarse,
                    num_ants: context.metafits_context.num_ants as u32,
                    num_timesteps: context.num_timesteps as u32,
                    num_pols: 1,
                    software: format!("Birli-{}", crate_version!()),
                    num_rows: (context.num_common_timesteps
                        * context.metafits_context.num_baselines)
                        as u32,
                    aoflagger_version: aoflagger_version.clone(),
                    aoflagger_strategy: aoflagger_strategy.clone(),
                };
                FlagFileSet::write_primary_hdu(&mut fptr, &header, Some(gpubox_id)).unwrap();
            }
        }

        for (&gpubox_id, path) in &gpubox_paths {
            let mut flag_fptr = FitsFile::open(path).unwrap();
            let hdu = flag_fptr.primary_hdu().unwrap();

            let obsid: i32 = get_required_fits_key!(&mut flag_fptr, &hdu, "OBSID").unwrap();
            assert_eq!(obsid, context.metafits_context.obs_id as i32);

            let header_gps_start: String =
                get_required_fits_key!(&mut flag_fptr, &hdu, "GPSSTART").unwrap();
            assert_abs_diff_eq!(header_gps_start.parse::<f64>().unwrap(), gps_start);

            let num_chans: i32 = get_required_fits_key!(&mut flag_fptr, &hdu, "NCHANS").unwrap();
            assert_eq!(
                num_chans,
                context.metafits_context.num_corr_fine_chans_per_coarse as i32
            );

            let num_ants: i32 = get_required_fits_key!(&mut flag_fptr, &hdu, "NANTENNA").unwrap();
            assert_eq!(num_ants, context.metafits_context.num_ants as i32);

            let num_scans: i32 = get_required_fits_key!(&mut flag_fptr, &hdu, "NSCANS").unwrap();
            assert_eq!(num_scans, context.num_timesteps as i32);

            let gpubox_no: i32 = get_required_fits_key!(&mut flag_fptr, &hdu, "GPUBOXNO").unwrap();
            assert_eq!(gpubox_no, gpubox_id as i32);
        }
    }

    #[test]
    fn test_write_occupancy_legacy() {
        let test_dir = Path::new("tests/data/1247842824_flags/");

        let corr_ctx = CorrelatorContext::new(
            test_dir.join("1247842824.metafits"),
            &[test_dir.join("1247842824_20190722150008_gpubox01_00.fits")],
        )
        .unwrap();
        let meta_ctx = &corr_ctx.metafits_context;
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        // TODO: only select a subset of baselines
        #[rustfmt::skip]
        let ant_pairs = vec![
            (0, 0), (0, 1), (0, 2), (0, 3),
            (1, 1), (1, 2), (1, 3),
            (2, 2), (2, 3),
            (3, 3),
        ];
        vis_sel.baseline_idxs = meta_ctx
            .baselines
            .iter()
            .enumerate()
            .filter(|(_, bl)| ant_pairs.contains(&(bl.ant1_index, bl.ant2_index)))
            .map(|(bl_idx, _)| bl_idx)
            .collect();
        vis_sel.timestep_range = 0..5;
        vis_sel.coarse_chan_range = 0..3;
        let fine_chans_per_coarse = meta_ctx.num_corr_fine_chans_per_coarse;

        let num_coarse_chans = vis_sel.coarse_chan_range.len();
        let num_timesteps = vis_sel.timestep_range.len();
        // let (num_timesteps, num_chans, num_baselines) = vis_sel.get_shape(fine_chans_per_coarse);
        let ant_indices = ant_indices(&vis_sel.get_ant_pairs(meta_ctx));
        let num_ants = ant_indices.len();

        let tmp_dir = tempdir().unwrap();
        let flag_template = tmp_dir.path().join("Flagfile%%.mwaf");
        let flag_template = flag_template.to_str().unwrap();

        let gpubox_ids: Vec<usize> = vis_sel
            .coarse_chan_range
            .clone()
            .map(|chan| corr_ctx.coarse_chans[chan].gpubox_number)
            .collect();

        #[cfg(feature = "aoflagger")]
        let (aoflagger_version, aoflagger_strategy) = {
            let mut major = 0;
            let mut minor = 0;
            let mut subminor = 0;
            unsafe {
                aoflagger_sys::cxx_aoflagger_new().GetVersion(
                    &mut major,
                    &mut minor,
                    &mut subminor,
                );
            }
            (
                Some(format!("v{major}.{minor}.{subminor}")),
                Some("amazing_flagging_strategy.lua".to_string()),
            )
        };
        #[cfg(not(feature = "aoflagger"))]
        let (aoflagger_version, aoflagger_strategy) = (None, None);

        let mut flag_set = FlagFileSet::new(
            flag_template,
            &corr_ctx,
            &vis_sel,
            aoflagger_version,
            aoflagger_strategy,
        )
        .unwrap();

        let mut flag_ctx = FlagContext::blank_from_dimensions(
            num_timesteps,
            num_coarse_chans,
            fine_chans_per_coarse,
            num_ants,
        );
        let int_time = 1e-3 * meta_ctx.corr_int_time_ms as f32;
        dbg!(int_time);
        flag_ctx.flag_init = int_time;
        flag_ctx.flag_end = int_time;
        flag_ctx.timestep_flags[1] = true;
        flag_ctx.coarse_chan_flags[1] = true;
        flag_ctx.fine_chan_flags[1] = true;
        flag_ctx.antenna_flags[1] = true;
        flag_ctx.flag_dc = true;
        flag_ctx.autos = true;

        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        flag_ctx
            .set_flags(
                flag_array.view_mut(),
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &ant_pairs,
            )
            .unwrap();

        flag_set.write_flag_array(flag_array.view(), false).unwrap();
        flag_set.finalise().unwrap();

        let flag_set =
            FlagFileSet::open(flag_template, &gpubox_ids, MWAVersion::CorrLegacy).unwrap();

        let (_, ch_occ) = flag_set.read_ch_occ().unwrap();
        for (cc_idx, occ) in ch_occ.chunks(fine_chans_per_coarse).enumerate() {
            if cc_idx == 1 {
                assert!(occ.iter().all(|&x| abs_diff_eq!(x, 1.)));
            } else {
                assert_abs_diff_eq!(occ[0], 0.76);
                assert_abs_diff_eq!(occ[1], 1.);
                assert_abs_diff_eq!(occ[2], 0.76);
            }
        }

        let (result_ant_pairs, _, bl_occ) = flag_set.read_bl_occ().unwrap();

        assert_eq!(
            result_ant_pairs,
            ant_pairs
                .iter()
                .map(|(a, b)| (*a as u32, *b as u32))
                .collect_vec()
        );
        for ((ant1, ant2), occ) in izip!(result_ant_pairs, bl_occ.axis_iter(Axis(1))) {
            if ant1 == ant2 || ant1 == 1 || ant2 == 1 {
                // assert_abs_diff_eq!(occ[0], 1.);
                // assert!(occ.)
                assert!(occ.iter().all(|&x| abs_diff_eq!(x, 1.)));
            } else {
                assert_abs_diff_eq!(occ[0], 0.2125);
                assert_abs_diff_eq!(occ[1], 1.);
                assert_abs_diff_eq!(occ[2], 0.2125);
            }
        }

        // TODO: finish test for bl_occ
    }

    #[test]
    fn test_read_flags_raw() {
        let test_dir = Path::new("tests/data/1247842824_flags/");

        let context = CorrelatorContext::new(
            test_dir.join("1247842824.metafits"),
            &[test_dir.join("1247842824_20190722150008_gpubox01_00.fits")],
        )
        .unwrap();
        let vis_sel = VisSelection::from_mwalib(&context).unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        // The old cotter mwaf files don't work, because the new format has new
        // keys.
        let filename_template = &test_dir.join("FlagfileCotterMWA%%.mwaf");
        let result = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        );
        assert!(result.is_err());

        // Generate new flags for testing.
        let temp_dir = tempdir().unwrap();
        let template = temp_dir
            .path()
            .join("FlagfileMWA%%.mwaf")
            .to_str()
            .unwrap()
            .to_string();
        let mut flag_file_set =
            FlagFileSet::new(&template, &context, &vis_sel, None, None).unwrap();
        let (num_timesteps, num_channels, num_baselines) = (2, 128, 8256);
        let mut flag_array = Array3::from_elem((num_timesteps, num_channels, num_baselines), false);
        flag_array.slice_mut(s![.., 0..8, ..]).fill(true);
        flag_array.slice_mut(s![.., 64, ..]).fill(true);
        flag_array.slice_mut(s![.., 120..128, ..]).fill(true);
        flag_array.slice_mut(s![0, 31, 4159]).fill(true);
        flag_file_set
            .write_flag_array(flag_array.view(), false)
            .unwrap();
        flag_file_set.finalise().unwrap();

        let flag_file_set = FlagFileSet::open(&template, &[1], MWAVersion::CorrLegacy).unwrap();
        // The shape of this array is (num_timesteps, num_baselines, num_channels)
        let disk_flags = flag_file_set.read_flags().unwrap();
        for (i_timestep, disk_flags) in disk_flags.outer_iter().enumerate() {
            for (i_baseline, disk_flags) in disk_flags.outer_iter().enumerate() {
                for (i_channel, disk_flag) in disk_flags.iter().enumerate() {
                    if [
                        0, 1, 2, 3, 4, 5, 6, 7, 64, 120, 121, 122, 123, 124, 125, 126, 127,
                    ]
                    .contains(&i_channel)
                    {
                        assert_eq!(
                            *disk_flag, 1, "timestep {i_timestep}, baseline {i_baseline}, {i_channel} channel was not flagged when it should've been!"
                        );
                        continue;
                    } else if i_timestep == 0 && i_channel == 31 && i_baseline == 4159 {
                        assert_eq!(*disk_flag, 1);
                    } else {
                        assert_eq!(*disk_flag, 0);
                    }
                }
            }
        }
    }
}
