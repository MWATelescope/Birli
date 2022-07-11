//! Items related to the reading and writing of the FITS-based MWA Flag file format.
//!
//! # MWAF Format
//!
//! Similar to the GPU Fits format, mwaf files come in a set for each observation, and there is one
//! .mwaf file per gpubox (coarse channel). This file contains a binary table of all the flags for
//! that coarse channel. There is one row for each timestep-baseline combination, and there is only
//! one column. Each cell in the table contains a binary vector of flags for each fine channel in
//! the coarse channel.

use std::{
    ops::Range,
    path::{Path, PathBuf},
};

use clap::crate_version;
use fitsio::{tables::ColumnDataType, tables::ColumnDescription, FitsFile};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::Itertools;
use marlu::{fitsio, fitsio_sys, mwalib, ndarray, rayon};
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
pub struct FlagFileHeader {
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
    /// The width of each fine channel mask vector in bytes, or the `NAXIS1` key from the table hdu
    pub bytes_per_row: u32,
    /// The number of rows (timesteps × baselines), and the `NAXIS2` key from the table hdu.
    pub num_rows: u32,
    /// The version of aoflagger used to generate these flags.
    pub aoflagger_version: Option<String>,
    /// The aoflagger strategy used to generate these flags.
    pub aoflagger_strategy: Option<String>,
}

/// A group of .mwaf Files for the same observation
pub struct FlagFileSet {
    /// A map associating gpubox numbers with filenames as well as the flag
    /// counts for that gpubox (used to determine occupancy which is written to
    /// the mwaf files after all flags are written).
    pub(crate) gpuboxes: Vec<(usize, PathBuf, Vec<u64>)>,
    /// Flag file header.
    pub(crate) header: FlagFileHeader,
    /// The number of baseline rows that have been written.
    row_count: u64,
    /// The number of baseline rows that are expected when writing is complete.
    expected_rows: u64,
}

impl FlagFileSet {
    fn get_gpubox_filenames(
        mwa_version: MWAVersion,
        filename_template: &str,
        gpubox_ids: &[usize],
    ) -> Result<Vec<(usize, PathBuf, Vec<u64>)>, IOError> {
        let num_percents = match mwa_version {
            MWAVersion::CorrOldLegacy | MWAVersion::CorrLegacy => 2,
            _ => 3,
        };
        let re_percents = Regex::new(format!("%{{{},}}+", num_percents).as_str()).unwrap();

        if !re_percents.is_match(filename_template) {
            return Err(InvalidFlagFilenameTemplate {
                source_file: file!(),
                source_line: line!(),
                filename_template: String::from(filename_template),
            });
        }

        Ok(gpubox_ids
            .iter()
            .map(|&gpubox_id| {
                (
                    gpubox_id,
                    PathBuf::from(
                        re_percents
                            .replace(
                                filename_template,
                                format!("{:0width$}", gpubox_id, width = num_percents),
                            )
                            .to_string(),
                    ),
                    vec![],
                )
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
    /// * `coarse_chan_range` is the... range of coarse channel indices.
    ///
    /// * `first_timestep` is the index corresponding to the first timestep used
    ///   in the [`mwalib::CorrelatorContext`].
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
        timestep_range: Range<usize>,
        coarse_chan_range: Range<usize>,
        aoflagger_version: Option<String>,
        aoflagger_strategy: Option<String>,
    ) -> Result<Self, IOError> {
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
        let num_rows = timestep_range.len() * corr_ctx.metafits_context.num_baselines;

        let header = FlagFileHeader {
            version: "2.0".to_string(),
            obs_id: corr_ctx.metafits_context.obs_id,
            gps_start,
            num_channels: num_fine_per_coarse,
            num_ants: corr_ctx.metafits_context.num_ants as u32,
            num_timesteps: timestep_range.len() as u32,
            num_pols: 1,
            // TODO: use something like https://github.com/rustyhorde/vergen
            software: format!("Birli-{}", crate_version!()),
            bytes_per_row: num_fine_per_coarse / 8 + u32::from(num_fine_per_coarse % 8 != 0),
            num_rows: num_rows as u32,
            aoflagger_version,
            aoflagger_strategy,
        };

        gpuboxes
            .par_iter_mut()
            .try_for_each(|(gpubox_id, filename, counts)| {
                // `counts` is currently 0 capacity; make it the right length.
                counts.resize(num_fine_per_coarse as usize, 0);

                if filename.exists() {
                    std::fs::remove_file(&filename)?;
                }

                match FitsFile::create(Path::new(filename)).open() {
                    Ok(mut fptr) => Self::write_primary_hdu(&mut fptr, &header, Some(*gpubox_id)),
                    Err(fits_error) => Err(FitsOpen {
                        fits_error,
                        fits_filename: filename.display().to_string(),
                        source_file: file!(),
                        source_line: line!(),
                    }),
                }
            })?;

        Ok(Self {
            gpuboxes,
            header,
            row_count: 0,
            expected_rows: num_rows as u64,
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
            bytes_per_row: _,
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
        hdu.write_key(fptr, "NCHANS", *num_channels as u32)?;
        hdu.write_key(fptr, "NANTENNA", *num_ants as u32)?;
        hdu.write_key(fptr, "NSCANS", *num_timesteps as u32)?;
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

        let multi_progress = MultiProgress::with_draw_target(if draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        });
        // progress bars detailing the write progress of each mwaf file.
        let progress_bars = self
            .gpuboxes
            .iter()
            .map(|(gpubox_id, _, _)| {
                multi_progress.add({
                    let pb = ProgressBar::new((num_timesteps * num_baselines) as _)
                        .with_style(
                            ProgressStyle::default_bar()
                                .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                                .progress_chars("=> "),
                        )
                        .with_position(0)
                        .with_message(format!("mwaf {:03}", gpubox_id));
                    pb.set_draw_delta(500); // Not setting this slows things down because the PB updates so quickly!
                    pb
                })
            })
            .collect::<Vec<_>>();

        // Can't use rayon::join or rayon::scope because one thread is dedicated
        // to the multi-progress bar, and if the system only has 1 CPU thread,
        // no work gets done, and the program deadlocks! crossbeam is overkill
        // but works well.
        let scoped_threads_result = crossbeam_utils::thread::scope(|s| {
            s.spawn(move |_| {
                multi_progress.join().unwrap();
            });

            let h = s.spawn(|_| {
                self.gpuboxes
                    .par_iter_mut()
                    .zip(
                        flag_array
                            .axis_chunks_iter(Axis(1), num_fine_chans_per_coarse)
                            .into_par_iter(),
                    )
                    .zip(progress_bars)
                    .try_for_each(
                        |(((_, filename, counts), flag_coarse_chan_view), channel_progress)| {
                            Self::write_flag_array_inner(
                                filename,
                                flag_coarse_chan_view,
                                &channel_progress,
                                num_fine_chans_per_coarse,
                                counts,
                            )
                        },
                    )
            });
            h.join()
        });
        match scoped_threads_result {
            // Propagate anything that didn't panic.
            Ok(Ok(r)) => r?,
            // A panic. This ideally only happens because a programmer made a
            // mistake, but it could happen in drastic situations (e.g. hardware
            // failure).
            Err(_) | Ok(Err(_)) => panic!(
                "A panic occurred; the message should be above. You may need to disable progress bars."
            ),
        };

        self.row_count += (flag_array.len_of(Axis(0)) * flag_array.len_of(Axis(2))) as u64;

        Ok(())
    }

    /// This fallible function is run in parallel from `write_flag_array`.
    fn write_flag_array_inner(
        filename: &Path,
        flag_coarse_chan_view: ArrayView3<bool>,
        channel_progress: &ProgressBar,
        num_fine_chans_per_coarse: usize,
        counts: &mut [u64],
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
            for flag_baseline_view in flag_timestep_view.axis_iter(Axis(1)) {
                flag_cell
                    .iter_mut()
                    .zip_eq(flag_baseline_view)
                    .zip_eq(counts.iter_mut())
                    .for_each(|((a, &b), count)| {
                        *a = i8::from(b);
                        if b {
                            *count += 1;
                        }
                    });

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
    /// there are problems in creating the OCCUPANCY table in the mwaf files.
    ///
    pub fn finalise(self) -> Result<(), IOError> {
        // Complain if not enough rows were written.
        if self.row_count != self.expected_rows {
            return Err(IOError::MwafIncorrectFlagCount {
                count: self.row_count,
                expected: self.expected_rows,
            });
        }

        // Write the occupancy stats into a new HDU ("OCCUPANCY").
        self.gpuboxes
            .into_par_iter()
            .try_for_each(|(_, filename, counts)| {
                Self::finalise_inner(filename, counts, self.expected_rows as f64)
            })?;

        Ok(())
    }

    /// This fallible function is run in parallel from `finalise`.
    fn finalise_inner(
        filename: PathBuf,
        counts: Vec<u64>,
        total_row_count: f64,
    ) -> Result<(), IOError> {
        let mut fptr = FitsFile::edit(filename)?;
        let index_col = ColumnDescription::new("Index")
            .with_type(ColumnDataType::Int)
            .create()?;
        let count_col = ColumnDescription::new("Count")
            .with_type(ColumnDataType::Long)
            .create()?;
        let occ_col = ColumnDescription::new("Occupancy")
            .with_type(ColumnDataType::Double)
            .create()?;
        let hdu = fptr.create_table("OCCUPANCY", &[index_col, count_col, occ_col])?;
        hdu.write_col(
            &mut fptr,
            "Index",
            &(0..counts.len())
                .into_iter()
                .map(|i| i as u32)
                .collect::<Vec<u32>>(),
        )?;
        hdu.write_col(&mut fptr, "Count", &counts)?;
        hdu.write_col(
            &mut fptr,
            "Occupancy",
            &counts
                .into_iter()
                .map(|c| c as f64 / total_row_count)
                .collect::<Vec<f64>>(),
        )?;

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
        let bytes_per_row = get_required_fits_key!(fptr, &hdu1, "NAXIS1")?;
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
            bytes_per_row,
            num_rows,
            aoflagger_version,
            aoflagger_strategy,
        };
        let baselines = (header.num_ants * (header.num_ants + 1)) / 2;
        if header.num_rows != header.num_timesteps * baselines {
            return Err(ReadMwafError::Generic(format!(
                "File {}: Expected NSCANS * NANTENNA * (NANTENNA+1) / 2 = NAXIS2, found {} * {} != {}",
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
        for (_, fits_filename, _) in &gpuboxes {
            match FitsFile::open(Path::new(&fits_filename)) {
                Ok(mut fptr) => {
                    if header.is_none() {
                        header = Some(Self::read_header(&mut fptr).map_err(|e| e.to_string())?.0);
                    }
                }
                Err(fits_error) => {
                    return Err(FitsOpen {
                        fits_error,
                        fits_filename: fits_filename.display().to_string(),
                        source_file: file!(),
                        source_line: line!(),
                    }
                    .to_string())
                }
            }
        }

        Ok(Self {
            gpuboxes,
            header: header.unwrap(),
            row_count: 0,
            expected_rows: 0,
        })
    }

    #[cfg(test)]
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

        for (gpubox_id, fits_filename, _) in &gpuboxes {
            match FitsFile::open(Path::new(&fits_filename)) {
                Ok(mut fptr) => {
                    let hdu0 = fits_open_hdu!(&mut fptr, 0).unwrap();
                    let version = get_required_fits_key!(&mut fptr, &hdu0, "VERSION").unwrap();
                    let obs_id = get_required_fits_key!(&mut fptr, &hdu0, "GPSTIME").unwrap();
                    let num_channels = get_required_fits_key!(&mut fptr, &hdu0, "NCHANS").unwrap();
                    let num_ants = get_required_fits_key!(&mut fptr, &hdu0, "NANTENNA").unwrap();
                    let num_timesteps = get_required_fits_key!(&mut fptr, &hdu0, "NSCANS").unwrap();
                    let num_pols = get_required_fits_key!(&mut fptr, &hdu0, "NPOLS").unwrap();
                    let fgpubox_id: u32 =
                        get_required_fits_key!(&mut fptr, &hdu0, "GPUBOXNO").unwrap();
                    let software = get_required_fits_key!(&mut fptr, &hdu0, "COTVER").unwrap();
                    let fdate = get_required_fits_key!(&mut fptr, &hdu0, "COTVDATE").unwrap();

                    let hdu1 = fits_open_hdu!(&mut fptr, 1).unwrap();
                    let bytes_per_row = get_required_fits_key!(&mut fptr, &hdu1, "NAXIS1").unwrap();
                    let num_rows = get_required_fits_key!(&mut fptr, &hdu1, "NAXIS2").unwrap();

                    assert_eq!(*gpubox_id as u32, fgpubox_id);
                    gpubox_ids.push(fgpubox_id);

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
                            bytes_per_row,
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
                        fits_filename: fits_filename.display().to_string(),
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
            },
            date.unwrap(),
        ))
    }

    #[cfg(test)]
    pub(crate) fn read_flags(&self) -> Result<Array3<i8>, IOError> {
        let (_, filename, _) = &self.gpuboxes[0];
        let mut fptr = FitsFile::open(filename)?;
        let hdu = fits_open_hdu!(&mut fptr, 0)?;
        let num_timesteps = get_required_fits_key!(&mut fptr, &hdu, "NSCANS")?;
        let num_channels_per_mwaf: usize = get_required_fits_key!(&mut fptr, &hdu, "NCHANS")?;
        let total_num_channels = num_channels_per_mwaf * self.gpuboxes.len();
        let num_ants: usize = get_required_fits_key!(&mut fptr, &hdu, "NANTENNA")?;
        let num_baselines = (num_ants * (num_ants + 1)) / 2;
        let hdu = fits_open_hdu!(&mut fptr, 1)?;
        let num_rows: usize = get_required_fits_key!(&mut fptr, &hdu, "NAXIS2")?;
        assert_eq!(num_rows, num_timesteps * num_baselines);

        let mut out = Array3::zeros((num_timesteps, num_baselines, total_num_channels));
        drop(fptr);
        drop(hdu);

        let mut row_flags = Array1::zeros(num_channels_per_mwaf);
        for (i_gpubox, (_, filename, _)) in self.gpuboxes.iter().enumerate() {
            let mut fptr = FitsFile::open(filename)?;
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
                    fits_filename: String::from(&fptr.filename),
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
    };
    use approx::assert_abs_diff_eq;
    use fitsio::FitsFile;
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
                let tmp_dir = tempdir().unwrap();
                assert!(matches!(
                    FlagFileSet::new(
                        tmp_dir.path().join($template_suffix).to_str().unwrap(),
                        $corr_ctx,
                        $timestep_range,
                        $coarse_chan_range,
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
        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        let ok_gpuboxes = 0..2;
        let colliding_gpuboxes = 1..3;

        for gpubox_id in colliding_gpuboxes.clone() {
            let colliding_filename = tmp_dir
                .path()
                .join(format!("Flagfile{:03}.mwaf", gpubox_id));
            File::create(colliding_filename.to_str().unwrap()).unwrap();
        }

        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &context,
                0..1,
                ok_gpuboxes,
                None,
                None
            ),
            Ok(FlagFileSet { .. })
        ));
        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &context,
                0..1,
                colliding_gpuboxes,
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

        for (_, filename, _) in gpuboxes {
            let mut fptr = FitsFile::open(filename).unwrap();
            let hdu0 = fptr.primary_hdu().unwrap();

            let version: String = get_required_fits_key!(&mut fptr, &hdu0, "VERSION").unwrap();
            let obs_id: Option<i32> = get_optional_fits_key!(&mut fptr, &hdu0, "OBSID").unwrap();
            let gps_start: Option<f64> =
                get_optional_fits_key!(&mut fptr, &hdu0, "GPSSTART").unwrap();
            let num_channels: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NCHANS").unwrap();
            let num_ants: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NANTENNA").unwrap();
            let num_timesteps: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NSCANS").unwrap();
            let num_pols: i32 = get_required_fits_key!(&mut fptr, &hdu0, "NPOLS").unwrap();
            let gpubox_id = get_optional_fits_key!(&mut fptr, &hdu0, "GPUBOXNO").unwrap();
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
            assert_eq!(gpubox_id, gpubox_id.map(|g| g as u32));
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
                tmp_dir
                    .path()
                    .join(format!("Flagfile{:03}.mwaf", gpubox_id)),
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
                    bytes_per_row: num_fine_per_coarse / 8
                        + u32::from(num_fine_per_coarse % 8 != 0),
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
    fn test_read_flags_raw() {
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
        let timestep_range = *context.common_timestep_indices.first().unwrap()
            ..*context.common_timestep_indices.last().unwrap() + 1;
        let temp_dir = tempdir().unwrap();
        let template = temp_dir
            .path()
            .join("FlagfileMWA%%.mwaf")
            .to_str()
            .unwrap()
            .to_string();
        let mut flag_file_set =
            FlagFileSet::new(&template, &context, timestep_range, 0..1, None, None).unwrap();
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
