//! IO for the AOCAL format, a binary calibration solutions format used by
//! Andre Offringa's Calibrate software.

pub(crate) use super::error::ReadSolutionsError;

use std::{fs::File, io::BufReader, path::Path};

use crate::{
    marlu::{
        hifitime::{Duration, Epoch},
        num_complex::Complex,
        Jones,
    },
    ndarray::{prelude::*, Array3},
};
use byteorder::{LittleEndian, ReadBytesExt};
use marlu::hifitime;

/// All of the relevant information contained within an MWAOCAL .bin file.
#[derive(Clone)]
pub struct AOCalSols {
    /// a three dimensional array of jones matrix calibration solutions with
    /// dimensions `[timestep][tile][frequency]`
    pub di_jones: Array3<Jones<f64>>,

    /// The start timestamps of each timeblock used to produce these calibration
    /// solutions.
    pub start_timestamps: Vec<Epoch>,
    // pub obsid: Option<u32>,
}

impl AOCalSols {
    /// Reads an MWAOCAL .bin file and returns a struct of its' contents.
    ///
    /// # Errors
    ///
    /// Can throw [`ReadSolutionsError`] if the file format is not valid.
    pub fn read_andre_binary<T: AsRef<Path>>(file: T) -> Result<Self, ReadSolutionsError> {
        let file_str = file.as_ref().display().to_string();
        // open the file, wrapping the IO Error in one which displays the file path
        let mut bin_file = BufReader::new(File::open(file).map_err(|e| {
            std::io::Error::new(e.kind(), format!("{} when accessing {}", e, file_str))
        })?);
        // The first 7 bytes should be ASCII "MWAOCAL".
        let mwaocal_str = String::from_utf8(vec![
            bin_file.read_u8()?,
            bin_file.read_u8()?,
            bin_file.read_u8()?,
            bin_file.read_u8()?,
            bin_file.read_u8()?,
            bin_file.read_u8()?,
            bin_file.read_u8()?,
        ])
        .unwrap();
        if mwaocal_str.as_str() != "MWAOCAL" {
            return Err(ReadSolutionsError::AndreBinaryStr {
                file: file_str,
                got: mwaocal_str,
            });
        }
        for _ in 0..9 {
            match bin_file.read_u8()? {
                0 => (),
                v => {
                    return Err(ReadSolutionsError::AndreBinaryVal {
                        file: file_str,
                        expected: "0",
                        got: v.to_string(),
                    })
                }
            }
        }
        let num_timeblocks = bin_file.read_u32::<LittleEndian>()? as usize;
        let total_num_tiles = bin_file.read_u32::<LittleEndian>()? as usize;
        let total_num_fine_freq_chans = bin_file.read_u32::<LittleEndian>()? as usize;
        let num_polarisations = bin_file.read_u32::<LittleEndian>()? as usize;
        let t = bin_file.read_f64::<LittleEndian>()?;
        let start_time = if t.abs() < f64::EPSILON {
            None
        } else {
            Some(Epoch::from_gpst_seconds(t))
        };
        let t = bin_file.read_f64::<LittleEndian>()?;
        let end_time = if t.abs() < f64::EPSILON {
            None
        } else {
            Some(Epoch::from_gpst_seconds(t))
        };
        let mut di_jones_vec = vec![
            0.0;
            num_timeblocks
                * total_num_tiles
                * total_num_fine_freq_chans
                // Real and imag for each polarisation.
                * 2 * num_polarisations
        ];
        bin_file.read_f64_into::<LittleEndian>(&mut di_jones_vec)?;
        let di_jones_a4 = Array4::from_shape_vec(
            (
                num_timeblocks,
                total_num_tiles,
                total_num_fine_freq_chans,
                2 * num_polarisations,
            ),
            di_jones_vec,
        )
        .unwrap();
        let di_jones = di_jones_a4.map_axis(Axis(3), |view| {
            Jones::from([
                Complex::new(view[0], view[1]),
                Complex::new(view[2], view[3]),
                Complex::new(view[4], view[5]),
                Complex::new(view[6], view[7]),
            ])
        });

        Ok(Self {
            di_jones,
            // We'd really like to have the *actual* start times of each
            // timeblock, but this isn't possible with this format. Here is the
            // best effort.
            start_timestamps: match (start_time, end_time) {
                (None, None) => vec![],
                (Some(t), None) => vec![t],
                (Some(s), Some(e)) if s == e => vec![s],
                (Some(s), Some(e)) if s != e => {
                    let duration = e - s;
                    let average_duration_per_timeblock =
                        duration.in_seconds() / (num_timeblocks - 1) as f64;
                    let mut start_timestamps = vec![];
                    for i_timeblock in 0..num_timeblocks {
                        let start: Epoch = s + i_timeblock as f64
                            * Duration::from_f64(
                                average_duration_per_timeblock,
                                hifitime::Unit::Second,
                            );
                        start_timestamps.push(start);
                    }
                    start_timestamps
                }
                _ => {
                    panic!("start_time is None, but end_time is not. Something went wrong.")
                }
            },
            // obsid: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::types::TestJones;

    use super::*;

    #[test]
    fn test_read_invalid_andre_binary() {
        let file = "tests/data/1254670392_avg/1254670392.metafits";
        assert!(matches!(
            AOCalSols::read_andre_binary(file),
            Err(ReadSolutionsError::AndreBinaryStr { .. })
        ));
    }

    #[test]
    fn test_read_valid_andre_binary() {
        let file = "tests/data/1254670392_avg/1254690096.bin";
        let sols = AOCalSols::read_andre_binary(file).unwrap();
        assert_eq!(sols.di_jones.dim(), (1, 128, 768));

        // it's not really useful to test much past here at this point, but this is useful for debugging
        // for (chan_idx, chan_view) in sols.di_jones.axis_iter(Axis(2)).enumerate() {
        //     eprintln!("{:3}: {:?} ({:?})", chan_idx, chan_view[(0, 0)], chan_view[(0, 0)].norm_sqr());
        // }

        // values from cotter debugger at breakpoint https://github.com/MWATelescope/cotter/blob/0f99a09cb21721666d2ba328ab2257484b4cd183/applysolutionswriter.cpp#L33

        assert_abs_diff_eq!(
            TestJones::from(sols.di_jones[(0, 0, 0)]),
            // -exec p this._solutions[0 * this._nSolutionChannels + 0]
            TestJones::from([
                Complex::new(-0.05711880819681107, 0.8909723224701427),
                Complex::new(0., 0.),
                Complex::new(0., 0.),
                Complex::new(-0.3190681285208096, 0.8975262420831494),
            ])
        );
        assert_abs_diff_eq!(
            TestJones::from(sols.di_jones[(0, 127, 0)]),
            // -exec p this._solutions[127 * this._nSolutionChannels + 0]
            TestJones::from([
                Complex::new(-0.6573012593247391, -0.16069224007242128),
                Complex::new(0., 0.),
                Complex::new(0., 0.),
                Complex::new(-0.6117273561749572, -0.2292172271648513)
            ])
        );
        assert_abs_diff_eq!(
            TestJones::from(sols.di_jones[(0, 0, 767)]),
            // -exec p this._solutions[0 * this._nSolutionChannels + 767]
            TestJones::from([
                Complex::new(-0.8662694860974055, 0.9154347855961527),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-0.6907841560791172, 0.9504852607301532)
            ])
        );
        assert_abs_diff_eq!(
            TestJones::from(sols.di_jones[(0, 127, 767)]),
            // -exec p this._solutions[127 * this._nSolutionChannels + 767]
            TestJones::from([
                Complex::new(-0.985537320263424, -0.21269260191492487),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.5117965761905826, -0.689634655264256)
            ])
        );
    }
}
