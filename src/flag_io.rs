// use std::collections::BTreeMap;

use crate::error::BirliError;
// use cxx::UniquePtr;
// use fitsio::FitsFile;
use mwalib::CorrelatorContext;
use regex::Regex;

// A group of MWAF Files for the same observation
#[derive(Clone, Debug)]
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
// fptrs: BTreeMap<usize, FitsFile>,
}

impl FlagFileSet {
    pub fn new(
        context: &CorrelatorContext,
        filename_template: &str,
        gpubox_ids: Vec<usize>,
    ) -> Result<Self, BirliError> {
        let re_percents = Regex::new("%{2,}+").unwrap();

        if !re_percents.is_match(filename_template) {
            return Err(BirliError::InvalidFlagFilenameTemplateError {
                source_file: file!(),
                source_line: line!(),
                filename_template: String::from(filename_template),
            });
        }

        dbg!(&context);
        dbg!(&gpubox_ids);

        unimplemented!();

        // for gpubox_id in gpubox_ids {
        //     let gpubox_id_str = gpubox_id.to_str();
        //     let filename = re_percents.replace()
        // }
    }

    // pub fn open(file_pattern: String) -> Result<Self, BirliError> {
    //     unimplemented!();
    // }

    // pub fn write_headers(){
    //     unimplemented!();
    // }

    // pub fn write_row(ant1_index: usize, ant2_index: usize, flags: const bool*) {
    //     unimplemented!();
    // }
}

#[cfg(test)]
mod tests {
    use super::FlagFileSet;
    use mwalib::CorrelatorContext;

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

    #[test]
    #[should_panic]
    fn test_flagfileset_enforces_percents_in_filename_template() {
        let context = get_mwax_context();
        let gpubox_ids = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();
        let filename_template = "bad_filename_template.mwaf";
        FlagFileSet::new(&context, filename_template, gpubox_ids).unwrap();
    }
}
