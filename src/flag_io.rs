use crate::error::BirliError;
use std::collections::BTreeMap;
use std::path::Path;
// use cxx::UniquePtr;
use fitsio::FitsFile;
use mwalib::{CorrelatorContext, CorrelatorVersion};
use regex::Regex;

// A group of MWAF Files for the same observation
// #[derive(Clone, Debug)]
pub struct FlagFileSet {
    // // GPSTIME
    // obs_id: u32,
    // // NCHANS - number of correlator fine channels per flag file
    // num_channels_per_file: usize,
    // // NANTENNA
    // num_ants: usize,
    // // NSCANS
    // num_timesteps: usize,
    // // each GPUBOXNO
    // gpubox_ids: Vec<usize>,
    #[allow(dead_code)]
    gpubox_fptrs: BTreeMap<usize, FitsFile>,
}

impl FlagFileSet {
    pub fn new(
        context: &CorrelatorContext,
        filename_template: &str,
        gpubox_ids: &Vec<usize>,
    ) -> Result<Self, BirliError> {
        let num_percents = match context.corr_version {
            CorrelatorVersion::Legacy | CorrelatorVersion::OldLegacy => 2,
            _ => 3,
        };
        let re_percents = Regex::new(format!("%{{{},}}+", num_percents).as_str()).unwrap();

        if !re_percents.is_match(filename_template) {
            return Err(BirliError::InvalidFlagFilenameTemplate {
                source_file: file!(),
                source_line: line!(),
                filename_template: String::from(filename_template),
            });
        }

        let mut gpubox_fptrs: BTreeMap<usize, FitsFile> = BTreeMap::new();
        for &gpubox_id in gpubox_ids.iter() {
            let filename = re_percents.replace(
                filename_template,
                format!("{:0width$}", gpubox_id, width = num_percents),
            );
            match FitsFile::create(Path::new(&filename.to_string())).open() {
                Ok(fptr) => {
                    gpubox_fptrs.insert(gpubox_id, fptr);
                }
                Err(fits_error) => {
                    return Err(BirliError::FitsIO {
                        fits_error: fits_error,
                        fits_filename: filename.to_string(),
                        source_file: file!(),
                        source_line: line!(),
                    })
                }
            }
        }

        FlagFileSet::write_headers(&mut gpubox_fptrs, &context);

        Ok(FlagFileSet {
            gpubox_fptrs: gpubox_fptrs,
        })

        // for gpubox_id in gpubox_ids {
        //     let gpubox_id_str = gpubox_id.to_str();
        //     let filename = re_percents.replace()
        // }
    }

    // pub fn open(file_pattern: String) -> Result<Self, BirliError> {
    //     unimplemented!();
    // }

    pub fn write_headers(
        gpubox_fptrs: &mut BTreeMap<usize, FitsFile>,
        context: &CorrelatorContext,
    ) {
        for (&gpubox_id, mut fptr) in gpubox_fptrs.into_iter() {
            let primary_hdu = fptr.primary_hdu().unwrap();
            primary_hdu
                .write_key(&mut fptr, "GPSTIME", context.metafits_context.obs_id as i32)
                .unwrap();
            primary_hdu
                .write_key(
                    &mut fptr,
                    "NCHANS",
                    context.metafits_context.num_corr_fine_chans_per_coarse as i32,
                )
                .unwrap();
            primary_hdu
                .write_key(
                    &mut fptr,
                    "NANTENNA",
                    context.metafits_context.num_ants as i32,
                )
                .unwrap();
            primary_hdu
                .write_key(&mut fptr, "NSCANS", context.num_timesteps as i32)
                .unwrap();
            primary_hdu
                .write_key(&mut fptr, "GPUBOXNO", gpubox_id as i32)
                .unwrap();
        }
    }

    // pub fn write_row(ant1_index: usize, ant2_index: usize, flags: const bool*) {
    //     unimplemented!();
    // }
}

#[cfg(test)]
mod tests {
    use super::FlagFileSet;
    use crate::error::BirliError;
    use fitsio::FitsFile;
    use mwalib::{CorrelatorContext, _get_optional_fits_key, get_optional_fits_key};
    use regex::Regex;
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
        let mwax_gpubox_ids = mwax_context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();
        let mwa_ord_context = get_mwa_ord_context();
        let mwa_ord_gpubox_ids = mwa_ord_context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();

        macro_rules! test_percent_enforcement {
            ($context:expr, $template_suffix:expr, $gpubox_ids:expr, $expected:pat) => {
                let tmp_dir = tempdir().unwrap();
                assert!(matches!(
                    FlagFileSet::new(
                        $context,
                        tmp_dir.path().join($template_suffix).to_str().unwrap(),
                        $gpubox_ids
                    ),
                    $expected
                ))
            };
        };
        test_percent_enforcement!(
            &mwax_context,
            "mwax_no_percents.mwaf",
            &mwax_gpubox_ids,
            Err(BirliError::InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            &mwa_ord_context,
            "mwa_ord_no_percents.mwaf",
            &mwa_ord_gpubox_ids,
            Err(BirliError::InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            &mwax_context,
            "mwax_insufficient_percents_2_%%.mwaf",
            &mwax_gpubox_ids,
            Err(BirliError::InvalidFlagFilenameTemplate { .. })
        );
        test_percent_enforcement!(
            &mwa_ord_context,
            "mwa_ord_sufficient_percents_2_%%.mwaf",
            &mwa_ord_gpubox_ids,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            &mwax_context,
            "mwax_sufficient_percents_3_%%%.mwaf",
            &mwax_gpubox_ids,
            Ok(FlagFileSet { .. })
        );
        test_percent_enforcement!(
            &mwa_ord_context,
            "mwa_ord_sufficient_percents_3_%%%.mwaf",
            &mwax_gpubox_ids,
            Ok(FlagFileSet { .. })
        );
    }

    #[test]
    fn test_flagfileset_fails_with_existing() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
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
            FlagFileSet::new(&context, filename_template.to_str().unwrap(), &ok_gpuboxes).unwrap(),
            FlagFileSet { .. }
        ));
        assert!(matches!(
            FlagFileSet::new(
                &context,
                filename_template.to_str().unwrap(),
                &colliding_gpuboxes
            )
            .err(),
            Some(BirliError::FitsIO { .. })
        ));
    }

    #[test]
    fn test_write_headers() {
        let context = get_mwax_context();
        let gpubox_ids: Vec<usize> = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();

        let tmp_dir = tempdir().unwrap();
        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");

        {
            let mut gpubox_fptrs: BTreeMap<usize, FitsFile> = BTreeMap::new();
            let re_percents = Regex::new("%{3,}+").unwrap();
            for &gpubox_id in gpubox_ids.iter() {
                let filename = re_percents.replace(
                    filename_template.to_str().unwrap(),
                    format!("{:03}", gpubox_id),
                );
                gpubox_fptrs.insert(
                    gpubox_id,
                    FitsFile::create(Path::new(&filename.to_string()))
                        .open()
                        .unwrap(),
                );
            }

            FlagFileSet::write_headers(&mut gpubox_fptrs, &context);
        }

        for &gpubox_id in gpubox_ids.iter() {
            let mut flag_fptr = FitsFile::open(
                tmp_dir
                    .path()
                    .join(format!("Flagfile{:03}.mwaf", gpubox_id)),
            )
            .unwrap();
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
}
