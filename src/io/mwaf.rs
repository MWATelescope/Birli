//! Items related to the reading and writing of the FITS-based MWA Flag file format.
//!
//! # MWAF Format
//!
//! Similar to the GPU Fits format, mwaf files come in a set for each observation, and there is one
//! .mwaf file per gpubox (coarse channel). This file contains a binary table of all the flags for
//! that coarse channel. There is one row for each timestep-baseline combination, and there is only
//! one column. Each cell in the table containsa binary vector of flags for each fine channel in
//! the coarse channel.

use super::error::{
    IOError,
    IOError::{FitsIO, FitsOpen, InvalidFlagFilenameTemplate, InvalidGpuBox, MwafInconsistent},
};
use crate::{
    error::BirliError,
    ndarray::{Array3, Axis},
};
use clap::crate_version;
use fitsio::{
    hdu::FitsHdu,
    tables::{ColumnDataDescription, ColumnDataType, ConcreteColumnDescription},
    FitsFile,
};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use marlu::{
    fitsio, fitsio_sys,
    mwalib::{
        CorrelatorContext, MWAVersion, _get_required_fits_key, _open_hdu, fits_open_hdu,
        get_required_fits_key,
    },
};
use regex::Regex;
use std::collections::BTreeMap;
use std::path::Path;

/// flag metadata which for a particular flag file in the set.
pub struct FlagFileHeaders {
    /// The `VERSION` key from the primary hdu
    // TODO: what is this actually used for?
    pub version: String,
    /// The `GPSTIME` key from the primary hdu
    pub obs_id: u32,
    /// The number of correlator fine channels per flag file, and the `NCHANS` key from the primary hdu.
    pub num_channels: usize,
    /// Total number of antennas (tiles) in the array, and the `NANTENNA` key from the primary hdu
    pub num_ants: usize,
    /// Number of timesteps in the observation, and the `NSCANS` key from the primary hdu
    pub num_timesteps: usize,
    /// The `NPOLS` key from the primary hdu
    pub num_pols: usize,
    /// The `GPUBOXNO` key from the primary hdu
    pub gpubox_id: usize,
    /// The `COTVER` key from the primary hdu
    pub cotter_version: String,
    /// The `COTVDATE` key from the primary hdu
    pub cotter_version_date: String,
    /// The width of each fine channel mask vector in bytes, or the `NAXIS1` key from the table hdu
    pub bytes_per_row: usize,
    /// The number of rows (timesteps × baselines), and the `NAXIS2` key from the table hdu.
    pub num_rows: usize,
    // TODO: is it useful to output aoflagger version and strategy?
}

impl FlagFileHeaders {
    /// Construct the [`FlagFileHeaders`] struct corresponding to the provided mwalib context and
    /// gpubox id.
    pub fn from_gpubox_context(gpubox_id: usize, context: &CorrelatorContext) -> Self {
        let num_fine_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        Self {
            version: "1.0".to_string(),
            obs_id: context.metafits_context.obs_id,
            num_channels: context.metafits_context.num_corr_fine_chans_per_coarse,
            num_ants: context.metafits_context.num_ants,
            num_timesteps: context.num_common_timesteps,
            num_pols: 1,
            gpubox_id,
            cotter_version: format!("Birli-{}", crate_version!()),
            // TODO: use something like https://github.com/rustyhorde/vergen
            cotter_version_date: "2021-04-14".to_string(),
            bytes_per_row: num_fine_per_coarse / 8 + usize::from(num_fine_per_coarse % 8 != 0),
            num_rows: context.num_common_timesteps * context.metafits_context.num_baselines,
        }
    }
}

/// A group of .mwaf Files for the same observation
pub struct FlagFileSet {
    gpubox_fptrs: BTreeMap<usize, FitsFile>,
}

impl FlagFileSet {
    fn get_gpubox_filenames(
        mwa_version: MWAVersion,
        filename_template: &str,
        gpubox_ids: &[usize],
    ) -> Result<BTreeMap<usize, String>, IOError> {
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

        let gpubox_filenames: BTreeMap<usize, String> = gpubox_ids
            .iter()
            .map(|&gpubox_id| {
                (
                    gpubox_id,
                    re_percents
                        .replace(
                            filename_template,
                            format!("{:0width$}", gpubox_id, width = num_percents),
                        )
                        .to_string(),
                )
            })
            .collect();

        Ok(gpubox_filenames)
    }

    /// Create a new set of flag files
    ///
    /// `filename_template` is a template string which is expanded to the list of flag files in the
    /// set, by replacing the percent (`%`) characters with each coarse channel's zero-prefixed
    /// GPU box ID. This is to maintain backwards compatibility with Cotter.
    ///
    /// For MWA Ord (legacy, pre-2021) correlator observations, the GPU box ID is the two digit
    /// correlator channel host number corresponding to `corr_chan_number` in
    /// [`marlu::mwalib::CoarseChannel`]
    ///
    /// For MWAX correlator observations, the GPU box ID is the three-digit received channel number
    /// corresponding to `rec_chan_number` in [`marlu::mwalib::CoarseChannel`].
    ///
    /// Be sure to specify the correct number of percent characters.
    ///
    /// # Errors
    ///
    /// Will error with [`IOError::FitsOpen`] if there are files already present at the paths
    /// specified in filename template.
    ///
    /// Will error with [`IOError::InvalidFlagFilenameTemplate`] if an invalid flag filename
    /// template is provided (wrong number of percents).
    pub fn new(
        filename_template: &str,
        gpubox_ids: &[usize],
        mwa_version: MWAVersion,
    ) -> Result<Self, IOError> {
        let mut gpubox_fptrs: BTreeMap<usize, FitsFile> = BTreeMap::new();
        let gpubox_filenames =
            Self::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)?;
        for (gpubox_id, fits_filename) in gpubox_filenames {
            match FitsFile::create(Path::new(&fits_filename)).open() {
                Ok(fptr) => {
                    gpubox_fptrs.insert(gpubox_id, fptr);
                }
                Err(fits_error) => {
                    return Err(FitsOpen {
                        fits_error,
                        fits_filename,
                        source_file: file!(),
                        source_line: line!(),
                    })
                }
            }
        }

        Ok(Self { gpubox_fptrs })
    }

    /// Open an existing set of flag files, given an observation's MWA Version, the flag filename
    /// template, and a list of gpubox ids.
    ///
    /// # Errors
    ///
    /// will error with [`BirliError::FitsOpen`] if files don't exist
    ///
    pub fn open(
        filename_template: &str,
        gpubox_ids: &[usize],
        mwa_version: MWAVersion,
    ) -> Result<Self, IOError> {
        let mut gpubox_fptrs: BTreeMap<usize, FitsFile> = BTreeMap::new();
        let gpubox_filenames =
            Self::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)?;
        for (gpubox_id, fits_filename) in gpubox_filenames {
            match FitsFile::open(Path::new(&fits_filename)) {
                Ok(fptr) => {
                    gpubox_fptrs.insert(gpubox_id, fptr);
                }
                Err(fits_error) => {
                    return Err(FitsOpen {
                        fits_error,
                        fits_filename,
                        source_file: file!(),
                        source_line: line!(),
                    })
                }
            }
        }

        Ok(Self { gpubox_fptrs })
    }

    fn write_primary_hdu(
        fptr: &mut FitsFile,
        hdu: &FitsHdu,
        header: &FlagFileHeaders,
    ) -> Result<(), IOError> {
        hdu.write_key(fptr, "VERSION", header.version.to_string())?;
        hdu.write_key(fptr, "GPSTIME", header.obs_id)?;
        hdu.write_key(fptr, "NCHANS", header.num_channels as u32)?;
        hdu.write_key(fptr, "NANTENNA", header.num_ants as u32)?;
        hdu.write_key(fptr, "NSCANS", header.num_timesteps as u32)?;
        hdu.write_key(fptr, "NPOLS", header.num_pols as u32)?;
        hdu.write_key(fptr, "GPUBOXNO", header.gpubox_id as u32)?;
        hdu.write_key(fptr, "COTVER", header.cotter_version.to_string())?;
        hdu.write_key(fptr, "COTVDATE", header.cotter_version_date.to_string())?;
        Ok(())
    }

    fn write_table_hdu(
        fptr: &mut FitsFile,
        hdu: &FitsHdu,
        header: &FlagFileHeaders,
    ) -> Result<(), IOError> {
        hdu.write_key(fptr, "NAXIS1", header.bytes_per_row as u32)?;
        hdu.write_key(fptr, "NAXIS2", header.num_rows as u32)?;
        Ok(())
    }

    fn read_header(fptr: &mut FitsFile) -> Result<FlagFileHeaders, IOError> {
        let hdu0 = fits_open_hdu!(fptr, 0)?;
        let hdu1 = fits_open_hdu!(fptr, 1)?;
        let header = FlagFileHeaders {
            version: get_required_fits_key!(fptr, &hdu0, "VERSION")?,
            obs_id: get_required_fits_key!(fptr, &hdu0, "GPSTIME")?,
            num_channels: get_required_fits_key!(fptr, &hdu0, "NCHANS")?,
            num_ants: get_required_fits_key!(fptr, &hdu0, "NANTENNA")?,
            num_timesteps: get_required_fits_key!(fptr, &hdu0, "NSCANS")?,
            num_pols: get_required_fits_key!(fptr, &hdu0, "NPOLS")?,
            gpubox_id: get_required_fits_key!(fptr, &hdu0, "GPUBOXNO")?,
            cotter_version: get_required_fits_key!(fptr, &hdu0, "COTVER")?,
            cotter_version_date: get_required_fits_key!(fptr, &hdu0, "COTVDATE")?,
            bytes_per_row: get_required_fits_key!(fptr, &hdu1, "NAXIS1")?,
            num_rows: get_required_fits_key!(fptr, &hdu1, "NAXIS2")?,
        };
        let baselines = header.num_ants * (header.num_ants + 1) / 2;
        if header.num_rows != header.num_timesteps * baselines {
            return Err(MwafInconsistent {
                file: String::from(&fptr.filename),
                expected: "NSCANS * NANTENNA * (NANTENNA+1) / 2 = NAXIS2".to_string(),
                found: format!(
                    "{} * {} != {}",
                    header.num_timesteps, baselines, header.num_rows
                ),
            });
        }
        Ok(header)
    }

    /// Write flags to disk, given an observation's [`marlu::mwalib::CorrelatorContext`], and an ndarray
    /// of boolean flags for the observation into a file for each `gpubox_id`.
    ///
    /// The filename template should contain two or 3 percentage (`%`) characters which will be replaced
    /// by the gpubox id or channel number (depending on correlator type). See [`FlagFileSet::new`]
    ///
    /// # Errors
    ///
    /// Will error if the gpubox ids this flagset was initialized with is not contained in the
    /// provided [`marlu::mwalib::CorrelatorContext`].
    ///
    pub fn write_flag_array(
        &mut self,
        context: &CorrelatorContext,
        flag_array: &Array3<bool>,
        gpubox_ids: &[usize],
    ) -> Result<(), IOError> {
        let flag_dims = flag_array.dim();
        let num_timesteps = flag_dims.0;
        let num_baselines = flag_dims.2;
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        // a parent progress bar containing the progress bars associated with this function
        let multi_progress = MultiProgress::with_draw_target(ProgressDrawTarget::stderr());
        // a vector of progress bars for the visibility reading progress of each channel.
        let write_progress: Vec<ProgressBar> = gpubox_ids
            .iter()
            .map(|gpubox_id| {
                let channel_progress = multi_progress.add(
                    ProgressBar::new((num_timesteps * num_baselines) as _)
                        .with_style(
                            ProgressStyle::default_bar()
                                .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                                .progress_chars("=> "),
                        )
                        .with_message(format!("gpubox {:03}", gpubox_id)),
                );
                channel_progress.set_position(0);
                channel_progress
            })
            .collect();

        // The total progress.
        let total_progress = multi_progress.add(
            ProgressBar::new((num_timesteps * num_baselines * gpubox_ids.len()) as _)
                .with_style(
                    ProgressStyle::default_bar()
                        .template(
                            "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                        )
                        .progress_chars("=> "),
                )
                .with_position(0)
                .with_message("write flags"),
        );

        // TODO: parallelize
        // let thread_write_info: Vec<(usize, &mut FitsFile, ArrayBase<ViewRepr<&bool>, Dim<[usize; 3]>>)> = Vec::with_capacity(gpubox_ids.len());
        for (&gpubox_id, flag_coarse_chan_view, channel_progress) in izip!(
            gpubox_ids,
            flag_array.axis_chunks_iter(Axis(1), fine_chans_per_coarse),
            write_progress
        ) {
            match self.gpubox_fptrs.get_mut(&gpubox_id) {
                None => {
                    return Err(InvalidGpuBox {
                        expected: format!("{:?}", self.gpubox_fptrs.keys()),
                        found: format!("gpubox id {}", gpubox_id),
                    })
                }
                Some(fptr) => {
                    // thread_write_info.push((gpubox_id, fptr, flag_coarse_chan_view));

                    let primary_hdu = fits_open_hdu!(fptr, 0)?;
                    let header = FlagFileHeaders::from_gpubox_context(gpubox_id, context);
                    let fine_chans_per_coarse = header.num_channels;
                    Self::write_primary_hdu(fptr, &primary_hdu, &header)?;

                    let flags_colname = "FLAGS";
                    let table_hdu = fptr.create_table(
                        "EXTNAME".to_string(),
                        &[ConcreteColumnDescription {
                            name: flags_colname.to_string(),
                            data_type: ColumnDataDescription::vector(
                                ColumnDataType::Bit,
                                fine_chans_per_coarse,
                            ),
                        }],
                    )?;
                    Self::write_table_hdu(fptr, &table_hdu, &header)?;

                    let mut status = 0;
                    let mut row_idx = 0;

                    for flag_timestep_view in flag_coarse_chan_view.outer_iter() {
                        for flag_baseline_view in flag_timestep_view.axis_iter(Axis(1)) {
                            let mut flag_cell = flag_baseline_view
                                .iter()
                                .map(|&flag| i8::from(flag))
                                .collect::<Vec<_>>();
                            assert_eq!(flag_cell.len(), fine_chans_per_coarse);

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
                                fits_filename: String::from(&fptr.filename),
                                hdu_num: 1,
                                source_file: file!(),
                                source_line: line!(),
                            })?;

                            row_idx += 1;
                            channel_progress.inc(1);
                            total_progress.inc(1);
                        }
                    }
                }
            }
            channel_progress.finish();
        }
        total_progress.finish();

        Ok(())
    }
}

/// TODO: These are just for tests, and should be deprecated.
/// TODO: Why doesn't #[cfg(test)] work?
impl FlagFileSet {
    fn read_flags_raw(
        fptr: &mut FitsFile,
        flags_raw: &mut [i8],
        row_idx: Option<usize>,
    ) -> Result<(), BirliError> {
        let mut status = 0;
        let row_idx = row_idx.unwrap_or(0);
        unsafe {
            fitsio_sys::ffgcx(
                fptr.as_raw(),
                1,
                1 + row_idx as i64,
                1,
                flags_raw.len() as i64,
                flags_raw.as_mut_ptr(),
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

        Ok(())
    }

    /// Read raw flags and headers from disk, as a [`std::collections::BTreeMap`] mapping from each
    /// gpubox id to a tuple containing a [`FlagFileHeaders`] and the raw flags as a vector of bytes.
    ///
    /// # Errors
    ///
    /// Will error with [`BirliError::MwafInconsistent`] if the mwaf files have inconsistent headers
    ///
    pub fn read_chan_header_flags_raw(
        &mut self,
    ) -> Result<BTreeMap<usize, (FlagFileHeaders, Vec<i8>)>, BirliError> {
        let mut chan_header_flags_raw = BTreeMap::new();

        for (&gpubox_id, fptr) in &mut self.gpubox_fptrs {
            let header = Self::read_header(fptr)?;
            let num_baselines = header.num_ants * (header.num_ants + 1) / 2;
            let mut flags_raw: Vec<i8> =
                vec![0; header.num_timesteps * num_baselines * header.num_channels];
            for timestep_idx in 0..header.num_timesteps {
                for baseline_idx in 0..num_baselines {
                    let row_idx = (timestep_idx * num_baselines) + baseline_idx;
                    let start_bit_idx = row_idx * header.num_channels;
                    let end_bit_idx = start_bit_idx + header.num_channels;
                    Self::read_flags_raw(
                        fptr,
                        &mut flags_raw[start_bit_idx..end_bit_idx],
                        Some(row_idx),
                    )?;
                }
            }
            chan_header_flags_raw.insert(gpubox_id, (header, flags_raw));
        }

        Ok(chan_header_flags_raw)
    }
}

#[cfg(test)]
mod tests {
    use super::{FlagFileHeaders, FlagFileSet};
    use crate::{
        io::error::IOError::{FitsOpen, InvalidFlagFilenameTemplate},
        test_common::{get_mwa_ord_context, get_mwax_context},
        VisSelection,
    };
    use fitsio::FitsFile;
    use marlu::{
        fitsio,
        mwalib::{
            CorrelatorContext, _get_optional_fits_key, _open_hdu, fits_open_hdu,
            get_optional_fits_key, MWAVersion,
        },
    };
    use ndarray::Axis;
    use std::collections::BTreeMap;
    use std::fs::File;
    use std::path::Path;
    use tempfile::tempdir;

    #[test]
    fn test_flagfileset_enforces_percents_in_filename_template() {
        let mwax_context = get_mwax_context();
        let mwax_gpubox_ids: Vec<_> = mwax_context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| mwax_context.coarse_chans[chan].gpubox_number)
            .collect();
        let mwa_ord_context = get_mwa_ord_context();
        let mwa_ord_gpubox_ids: Vec<_> = mwa_ord_context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| mwa_ord_context.coarse_chans[chan].gpubox_number)
            .collect();

        macro_rules! test_percent_enforcement {
            ($template_suffix:expr, $gpubox_ids:expr, $mwa_version:expr, $expected:pat) => {
                let tmp_dir = tempdir().unwrap();
                assert!(matches!(
                    FlagFileSet::new(
                        tmp_dir.path().join($template_suffix).to_str().unwrap(),
                        $gpubox_ids,
                        $mwa_version
                    ),
                    $expected
                ))
            };
        }
        test_percent_enforcement!(
            "mwax_no_percents.mwaf",
            &mwax_gpubox_ids,
            MWAVersion::CorrMWAXv2,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_no_percents.mwaf",
            &mwa_ord_gpubox_ids,
            MWAVersion::CorrLegacy,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwax_insufficient_percents_2_%%.mwaf",
            &mwax_gpubox_ids,
            MWAVersion::CorrMWAXv2,
            Err(InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_sufficient_percents_2_%%.mwaf",
            &mwa_ord_gpubox_ids,
            MWAVersion::CorrLegacy,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            "mwax_sufficient_percents_3_%%%.mwaf",
            &mwax_gpubox_ids,
            MWAVersion::CorrMWAXv2,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            "mwa_ord_sufficient_percents_3_%%%.mwaf",
            &mwax_gpubox_ids,
            MWAVersion::CorrLegacy,
            Ok(FlagFileSet { .. })
        );
    }

    #[test]
    fn test_flagfileset_new_fails_with_existing() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        let ok_gpuboxes = gpubox_ids[..1].to_vec();
        let colliding_gpuboxes = gpubox_ids[1..].to_vec();

        for gpubox_id in &colliding_gpuboxes {
            let colliding_filename = tmp_dir
                .path()
                .join(format!("Flagfile{:03}.mwaf", gpubox_id));
            File::create(colliding_filename.to_str().unwrap()).unwrap();
        }

        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &ok_gpuboxes,
                context.mwa_version
            ),
            Ok(FlagFileSet { .. })
        ));
        assert!(matches!(
            FlagFileSet::new(
                filename_template.to_str().unwrap(),
                &colliding_gpuboxes,
                context.mwa_version
            ),
            Err(FitsOpen { .. })
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

        assert!(matches!(
            FlagFileSet::open(
                filename_template.to_str().unwrap(),
                &gpuboxes,
                context.mwa_version
            ),
            Err(FitsOpen { .. })
        ));
    }

    #[test]
    fn test_read_headers() {
        let test_dir = Path::new("tests/data/1247842824_flags/");

        let context = CorrelatorContext::new(
            &test_dir.join("1247842824.metafits"),
            &[test_dir.join("1247842824_20190722150008_gpubox01_00.fits")],
        )
        .unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let filename_template = &test_dir.join("FlagfileCotterMWA%%.mwaf");
        let mut flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        for (&gpubox_id, fptr) in &mut flag_file_set.gpubox_fptrs {
            let header = FlagFileSet::read_header(fptr).unwrap();
            assert_eq!(header.obs_id, 1247842824);
            assert_eq!(header.num_channels, 128);
            assert_eq!(header.num_ants, 128);
            assert_eq!(header.num_timesteps, 2);
            assert_eq!(header.gpubox_id, gpubox_id);
            assert_eq!(header.cotter_version, "4.5");
            assert_eq!(header.cotter_version_date, "2020-08-10");
        }
    }

    #[test]
    fn test_write_primary_hdu() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

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

        {
            for (&gpubox_id, path) in &gpubox_paths {
                let mut fptr = FitsFile::create(path).open().unwrap();
                let primary_hdu = fits_open_hdu!(&mut fptr, 0).unwrap();
                FlagFileSet::write_primary_hdu(
                    &mut fptr,
                    &primary_hdu,
                    &FlagFileHeaders::from_gpubox_context(gpubox_id, &context),
                )
                .unwrap();
            }
        }

        for (&gpubox_id, path) in &gpubox_paths {
            let mut flag_fptr = FitsFile::open(path).unwrap();
            let hdu = flag_fptr.primary_hdu().unwrap();

            let gps_time: Option<i32> =
                get_optional_fits_key!(&mut flag_fptr, &hdu, "GPSTIME").unwrap();
            assert_eq!(gps_time.unwrap(), context.metafits_context.obs_id as i32);

            let num_chans: Option<i32> =
                get_optional_fits_key!(&mut flag_fptr, &hdu, "NCHANS").unwrap();
            assert_eq!(
                num_chans.unwrap(),
                context.metafits_context.num_corr_fine_chans_per_coarse as i32
            );

            let num_ants: Option<i32> =
                get_optional_fits_key!(&mut flag_fptr, &hdu, "NANTENNA").unwrap();
            assert_eq!(num_ants.unwrap(), context.metafits_context.num_ants as i32);

            let num_scans: Option<i32> =
                get_optional_fits_key!(&mut flag_fptr, &hdu, "NSCANS").unwrap();
            assert_eq!(num_scans.unwrap(), context.num_timesteps as i32);

            let gpubox_no: Option<i32> =
                get_optional_fits_key!(&mut flag_fptr, &hdu, "GPUBOXNO").unwrap();
            assert_eq!(gpubox_no.unwrap(), gpubox_id as i32);
        }
    }

    #[test]
    fn test_read_flags_raw() {
        let test_dir = Path::new("tests/data/1247842824_flags/");

        let context = CorrelatorContext::new(
            &test_dir.join("1247842824.metafits"),
            &[test_dir.join("1247842824_20190722150008_gpubox01_00.fits")],
        )
        .unwrap();

        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let filename_template = &test_dir.join("FlagfileCotterMWA%%.mwaf");
        let mut flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        let chan_hdr_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();

        assert_eq!(chan_hdr_flags_raw.keys().len(), 1);
        let (chan1_header, chan1_flags_raw) = &chan_hdr_flags_raw[&1];
        assert!(!chan1_flags_raw.is_empty());

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;

        let tests = [
            (0, 2, 103, i8::from(false)),
            (0, 2, 104, i8::from(true)),
            (0, 2, 105, i8::from(true)),
            (0, 2, 106, i8::from(false)),
            (1, 2, 103, i8::from(false)),
            (1, 2, 104, i8::from(true)),
            (1, 2, 105, i8::from(true)),
            (1, 2, 106, i8::from(false)),
        ];
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in &tests {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let row = &chan1_flags_raw[(row_idx * chan1_header.num_channels)
                ..((row_idx + 1) * chan1_header.num_channels)];
            assert_eq!(
                &row[*fine_chan_idx], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, row {:?}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, row
            );
        }
    }

    #[test]
    fn test_write_flag_array() {
        let corr_ctx = get_mwax_context();

        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        let gpubox_ids = corr_ctx.coarse_chans[vis_sel.coarse_chan_range.clone()]
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect::<Vec<_>>();

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();

        let mut idx = 0;

        for (coarse_chan_idx, mut flag_coarse_chan_view) in flag_array
            .axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse)
            .enumerate()
        {
            for (timestep_idx, mut flag_timestep_view) in
                flag_coarse_chan_view.outer_iter_mut().enumerate()
            {
                for (baseline_idx, mut flag_baseline_view) in
                    flag_timestep_view.axis_iter_mut(Axis(1)).enumerate()
                {
                    for (fine_chan_idx, fine_chan_flag) in flag_baseline_view.iter_mut().enumerate()
                    {
                        *fine_chan_flag = 1 << fine_chan_idx & idx != 0;
                        println!(
                            "{} (cc {}, ts {}, bl {} fc {}) = {}",
                            idx,
                            coarse_chan_idx,
                            timestep_idx,
                            baseline_idx,
                            fine_chan_idx,
                            fine_chan_flag
                        );
                    }
                    idx = (idx + 1) % (1 << fine_chans_per_coarse);
                }
            }
        }

        let mut flag_file_set = FlagFileSet::new(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();
        flag_file_set
            .write_flag_array(&corr_ctx, &flag_array, &gpubox_ids)
            .unwrap();
        drop(flag_file_set);

        let mut flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            corr_ctx.mwa_version,
        )
        .unwrap();

        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();

        assert_eq!(chan_header_flags_raw.keys().len(), 2);
        let (chan1_header, chan1_flags_raw) = chan_header_flags_raw.get(&117).unwrap();
        dbg!(chan1_flags_raw);

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, corr_ctx.num_common_timesteps);
        assert_eq!(num_baselines, corr_ctx.metafits_context.num_baselines);
        assert_eq!(chan1_header.num_channels, fine_chans_per_coarse);
        assert_eq!(
            chan1_flags_raw.len(),
            chan1_header.num_timesteps * num_baselines * chan1_header.num_channels
        );

        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in &[
            (0, 0, 0, i8::from(false)),
            (0, 0, 1, i8::from(false)),
            (0, 1, 0, i8::from(true)),
            (0, 1, 1, i8::from(false)),
            (0, 2, 0, i8::from(false)),
            (0, 2, 1, i8::from(true)),
            (1, 0, 0, i8::from(true)),
            (1, 0, 1, i8::from(true)),
            (1, 1, 0, i8::from(false)),
            (1, 1, 1, i8::from(false)),
            (1, 2, 0, i8::from(true)),
            (1, 2, 1, i8::from(false)),
        ] {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx * fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                &chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }
}
