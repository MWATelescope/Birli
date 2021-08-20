//! Items related to the reading and writing of the FITS-based MWA Flag file format.
//!
//! # MWAF Format
//!
//! Similar to the GPUFits format, mwaf files come in a set for each observation, and there is one
//! .mwaf file per gpubox (coarse channel). This file contains a binary table of all the flags for
//! that coarse channel. There is one row for each timestep-baseline combination, and there is only
//! one column. Each cell in the table containsa binary vector of flags for each fine channel in
//! the coarse channel.

use super::error::{
    IOError,
    IOError::{FitsIO, FitsOpen, InvalidFlagFilenameTemplate, InvalidGpuBox, MwafInconsistent},
};
use crate::cxx_aoflagger::ffi::CxxFlagMask;
use crate::error::BirliError;
use clap::crate_version;
use cxx::UniquePtr;
use fitsio::hdu::FitsHdu;
use fitsio::tables::{ColumnDataDescription, ColumnDataType, ConcreteColumnDescription};
use fitsio::FitsFile;
use mwa_rust_core::{fitsio, fitsio_sys, mwalib};
use mwalib::{
    CorrelatorContext, MWAVersion, _get_required_fits_key, _open_hdu, fits_open_hdu,
    get_required_fits_key,
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
        FlagFileHeaders {
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

// TODO: can this just take metafits context instead of a full context?
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
    /// GPUBox ID. This is to maintain backwards compatibility with Cotter.
    ///
    /// For MWA Ord (legacy, pre-2021) correlator observations, the GPUBox ID is the two digit
    /// correlator channel host number corresponding to [`mwalib::CoarseChannel.corr_chan_number`]
    ///
    /// For MWAX correlator observations, the GPUBox ID is the three-digit received channel number
    /// corresponding to [`mwalib::CoarseChannel.rec_chan_number`].
    ///
    /// Be sure to specify the correct number of percent characters.
    ///
    /// # Errors
    ///
    /// Will error with [`BirliError::FitsOpen`] if there are files already present at the paths
    /// specified in filename template.
    ///
    /// Will error with [`BirliError::InvalidFlagFilenameTemplate`] if an invalid flag filename
    /// template is provided (wrong number of percents).
    pub fn new(
        filename_template: &str,
        gpubox_ids: &[usize],
        mwa_version: MWAVersion,
    ) -> Result<Self, IOError> {
        let mut gpubox_fptrs: BTreeMap<usize, FitsFile> = BTreeMap::new();
        let gpubox_filenames =
            FlagFileSet::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)?;
        for (gpubox_id, fits_filename) in gpubox_filenames.into_iter() {
            match FitsFile::create(Path::new(&fits_filename.to_string())).open() {
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

        Ok(FlagFileSet { gpubox_fptrs })
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
            FlagFileSet::get_gpubox_filenames(mwa_version, filename_template, gpubox_ids)?;
        for (gpubox_id, fits_filename) in gpubox_filenames.into_iter() {
            match FitsFile::open(Path::new(&fits_filename.to_string())) {
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

        Ok(FlagFileSet { gpubox_fptrs })
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

    /// Write flags to disk, given an observation's [`mwalib::CorrelatorContext`], and a
    /// vector of [`CxxFlagMask`]s for each baseline in the observation.
    ///
    /// The filename template should contain two or 3 percentage (`%`) characters which will be replaced
    /// by the gpubox id or channel number (depending on correlator type). See [`FlagFileSet::new`]
    ///
    /// # Errors
    ///
    /// Will error if the gpubox ids this flagset was initialized with is not contained in the
    /// provided [`mwalib::CorrelatorContext`].
    ///
    pub fn write_baseline_flagmasks(
        &mut self,
        context: &CorrelatorContext,
        baseline_flagmasks: &[UniquePtr<CxxFlagMask>],
        img_coarse_chan_idxs: &[usize],
    ) -> Result<(), IOError> {
        // let gpubox_chan_numbers: BTreeMap<usize, usize> = context
        //     .coarse_chans
        //     .iter()
        //     .map(|chan| (chan.gpubox_number, chan.corr_chan_number))
        //     .collect();
        let gpubox_to_chan_idx: BTreeMap<usize, usize> = context
            .coarse_chans
            .iter()
            .enumerate()
            .map(|(chan_idx, chan)| (chan.gpubox_number, chan_idx))
            .collect();

        for (&gpubox_id, fptr) in self.gpubox_fptrs.iter_mut() {
            let primary_hdu = fits_open_hdu!(fptr, 0)?;
            let header = FlagFileHeaders::from_gpubox_context(gpubox_id, context);
            let num_baselines = header.num_ants * (header.num_ants + 1) / 2;
            let num_fine_chans_per_coarse = header.num_channels;
            FlagFileSet::write_primary_hdu(fptr, &primary_hdu, &header)?;
            // let chan_number = match gpubox_chan_numbers.get(&gpubox_id) {
            //     Some(chan_number) => chan_number,
            //     None => {
            //         return Err(InvalidGpuBox {
            //             expected: format!("{:?}", gpubox_chan_numbers.keys()),
            //             found: format!("{}", gpubox_id),
            //         })
            //     }
            // };
            let coarse_chan_idx = match gpubox_to_chan_idx.get(&gpubox_id) {
                Some(chan_number) => chan_number,
                None => {
                    return Err(InvalidGpuBox {
                        expected: format!("{:?}", gpubox_to_chan_idx.keys()),
                        found: format!("{}", gpubox_id),
                    })
                }
            };
            let img_coarse_chan_idx = match img_coarse_chan_idxs
                .iter()
                .position(|idx| idx == coarse_chan_idx)
            {
                Some(img_coarse_chan_idx) => img_coarse_chan_idx,
                None => {
                    return Err(InvalidGpuBox {
                        expected: format!("{:?}", img_coarse_chan_idxs),
                        found: format!("{}", coarse_chan_idx),
                    })
                }
            };
            let flags_colname = "FLAGS";
            let table_hdu = fptr.create_table(
                "EXTNAME".to_string(),
                &[ConcreteColumnDescription {
                    name: flags_colname.to_string(),
                    data_type: ColumnDataDescription::vector(
                        ColumnDataType::Bit,
                        num_fine_chans_per_coarse,
                    ),
                }],
            )?;
            FlagFileSet::write_table_hdu(fptr, &table_hdu, &header)?;

            let mut status = 0;
            for (baseline_idx, flagmask) in baseline_flagmasks.iter().enumerate() {
                assert!(baseline_idx < num_baselines);
                // Flag buffer can be thought of as a 2D buffer of bools,
                // indexed by channel (y), then timestep (x). flag stride is the width
                let flag_buffer = flagmask.Buffer();
                let flag_stride = flagmask.HorizontalStride();
                let num_timesteps = flagmask.Width();
                let num_channels = flagmask.Height();
                assert_eq!(
                    num_channels,
                    img_coarse_chan_idxs.len() * num_fine_chans_per_coarse
                );
                for img_timestep_idx in 0..num_timesteps {
                    let row_idx = (img_timestep_idx * num_baselines) + baseline_idx;
                    // All flags for img_timestep_idx
                    let timestep_flags = flag_buffer
                        .iter()
                        .skip(img_timestep_idx)
                        .step_by(flag_stride);

                    // All flags for timestep_idx and img_coarse_chan_idx
                    let mut cell: Vec<i8> = timestep_flags
                        .skip(img_coarse_chan_idx * num_fine_chans_per_coarse)
                        .take(num_fine_chans_per_coarse)
                        .map(|&flag| i8::from(flag))
                        .collect();
                    unsafe {
                        fitsio_sys::ffpclx(
                            fptr.as_raw(),
                            1,
                            1 + row_idx as i64,
                            1,
                            cell.len() as i64,
                            cell.as_mut_ptr(),
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
                }
            }
        }
        Ok(())
    }
}

/// TODO: These are just for tests, and should be deprecated.
/// TODO: Why doesn't #[cfg(test)] work?
impl FlagFileSet {
    // pub fn read_validated_header(
    //     context: &CorrelatorContext,
    //     fptr: &mut FitsFile,
    // ) -> Result<FlagFileHeaders, BirliError> {
    //     let headers = FlagFileSet::read_header(fptr)?;
    //     let header_baselines = headers.num_ants * (headers.num_ants + 1) / 2;
    //     if header_baselines != context.metafits_context.num_baselines {
    //         return Err(BirliError::MwafInconsistent {
    //             file: String::from(&fptr.filename),
    //             expected: "NANTENNA * (NANTENNA+1) / 2 = context.metafits_context.num_baselines"
    //                 .to_string(),
    //             found: format!(
    //                 "{} != {}",
    //                 header_baselines, context.metafits_context.num_baselines
    //             ),
    //         });
    //     };

    //     // TODO: check NSCANS?
    //     // if headers.num_timesteps > context.num_common_timesteps {
    //     //     return Err(BirliError::MwafInconsistent {
    //     //         file: String::from(&fptr.filename),
    //     //         expected: "NSCANS <= context.num_common_timesteps".to_string(),
    //     //         found: format!("{} > {}", headers.num_timesteps, context.num_common_timesteps),
    //     //     });
    //     // };

    //     if headers.bytes_per_row * 8 < context.metafits_context.num_corr_fine_chans_per_coarse {
    //         return Err(BirliError::MwafInconsistent {
    //             file: String::from(&fptr.filename),
    //             expected: "headers.bytes_per_row * 8 >= context.metafits_context.num_corr_fine_chans_per_coarse".to_string(),
    //             found: format!("{} < {}", headers.bytes_per_row, context.metafits_context.num_corr_fine_chans_per_coarse),
    //         });
    //     }

    //     Ok(headers)
    // }

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

        for (&gpubox_id, fptr) in self.gpubox_fptrs.iter_mut() {
            let header = FlagFileSet::read_header(fptr)?;
            let num_baselines = header.num_ants * (header.num_ants + 1) / 2;
            let mut flags_raw: Vec<i8> =
                vec![0; header.num_timesteps * num_baselines * header.num_channels];
            for timestep_idx in 0..header.num_timesteps {
                for baseline_idx in 0..num_baselines {
                    let row_idx = (timestep_idx * num_baselines) + baseline_idx;
                    let start_bit_idx = row_idx * header.num_channels;
                    let end_bit_idx = start_bit_idx + header.num_channels;
                    FlagFileSet::read_flags_raw(
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
    use crate::cxx_aoflagger::ffi::{cxx_aoflagger_new, CxxFlagMask};
    use crate::io::error::IOError::{FitsOpen, InvalidFlagFilenameTemplate};
    use cxx::UniquePtr;
    use fitsio::FitsFile;
    use mwa_rust_core::{fitsio, mwalib};
    use mwalib::{
        CorrelatorContext, _get_optional_fits_key, _open_hdu, fits_open_hdu, get_optional_fits_key,
        MWAVersion,
    };
    use std::collections::BTreeMap;
    use std::fs::File;
    use std::path::Path;
    use tempfile::tempdir;

    // TODO: deduplicate this from lib.rs
    fn get_mwax_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
        let gpufits_paths = vec![
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    fn get_mwa_ord_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

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

        for gpubox_id in colliding_gpuboxes.iter() {
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

        for (&gpubox_id, mut fptr) in flag_file_set.gpubox_fptrs.iter_mut() {
            let header = FlagFileSet::read_header(&mut fptr).unwrap();
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
        for &gpubox_id in gpubox_ids.iter() {
            gpubox_paths.insert(
                gpubox_id,
                tmp_dir
                    .path()
                    .join(format!("Flagfile{:03}.mwaf", gpubox_id)),
            );
        }

        {
            for (&gpubox_id, path) in gpubox_paths.iter() {
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

        for (&gpubox_id, path) in gpubox_paths.iter() {
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

        for (_, fptr) in flag_file_set.gpubox_fptrs.iter_mut() {
            let table_hdu = fptr.hdu(1).unwrap();
            dbg!(table_hdu);
        }
        let chan_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();

        assert_eq!(chan_flags_raw.keys().len(), 1);
        let (chan1_header, chan1_flags_raw) = chan_flags_raw.get(&1).unwrap();
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
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests.iter() {
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
    fn test_write_baseline_flagmasks() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
            .collect();

        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        let mut idx = 0;
        let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let num_img_coarse_chans = img_coarse_chan_idxs.len();

        let width = context.num_common_timesteps;
        let height = num_img_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut baseline_flagmasks: Vec<UniquePtr<CxxFlagMask>> = context
            .metafits_context
            .baselines
            .iter()
            .map(|_| unsafe { aoflagger.MakeFlagMask(width, height, false) })
            .collect();

        for img_coarse_chan_idx in 0..num_img_coarse_chans {
            for img_timestep_idx in 0..width {
                for flag_mask_ptr in baseline_flagmasks.iter_mut() {
                    let flag_stride = flag_mask_ptr.HorizontalStride();
                    let flag_buf = flag_mask_ptr.pin_mut().BufferMut();
                    for fine_chan_idx in 0..num_fine_chans_per_coarse {
                        let flag_offset_y =
                            img_coarse_chan_idx * num_fine_chans_per_coarse + fine_chan_idx;
                        let flag_idx = flag_offset_y * flag_stride + img_timestep_idx;
                        assert!(flag_idx < flag_stride * height);
                        dbg!(flag_idx, fine_chan_idx, idx, 1 << fine_chan_idx & idx);
                        flag_buf[flag_idx] = 1 << fine_chan_idx & idx != 0;
                    }
                    idx = (idx + 1) % (1 << num_fine_chans_per_coarse);
                }
            }
        }

        let mut flag_file_set = FlagFileSet::new(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();
        flag_file_set
            .write_baseline_flagmasks(&context, &baseline_flagmasks, img_coarse_chan_idxs)
            .unwrap();
        drop(flag_file_set);

        let mut flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &gpubox_ids,
            context.mwa_version,
        )
        .unwrap();

        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();

        assert_eq!(chan_header_flags_raw.keys().len(), 2);
        let (chan1_header, chan1_flags_raw) = chan_header_flags_raw.get(&117).unwrap();
        dbg!(chan1_flags_raw);

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, context.num_common_timesteps);
        assert_eq!(num_baselines, context.metafits_context.num_baselines);
        assert_eq!(chan1_header.num_channels, num_fine_chans_per_coarse);
        assert_eq!(
            chan1_flags_raw.len(),
            chan1_header.num_timesteps * num_baselines * chan1_header.num_channels
        );

        let tests = [
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
        ];
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests.iter() {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx * num_fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                &chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }
}
