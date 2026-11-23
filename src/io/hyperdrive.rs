//! IO for the hyperdrive fits calibration solutions format.
//!
//! See <https://mwatelescope.github.io/mwa_hyperdrive/defs/cal_sols_hyp.html>
//!
//! Only the parts of the format that Birli needs to apply solutions are read:
//!
//! - the `SOLUTIONS` image, with dimensions `[timeblock][tile][chanblock][8]`,
//!   where the last axis is the real and imaginary parts of the XX, XY, YX and
//!   YY elements of each Jones matrix. Flagged tiles and chanblocks are NaN.
//! - the `Freq` column of the `CHANBLOCKS` table, the centre frequency \[Hz\]
//!   of each chanblock.
//! - the `TileName` column of the `TILES` table.
//! - the `Start` column of the `TIMEBLOCKS` table (GPS seconds).
//! - the `OBSID`, `PFB`, `D_GAINS`, `CABLELEN` and `GEOMETRY` keys of the
//!   primary HDU. The last four are only written when hyperdrive calibrated
//!   from raw MWA files, and say which corrections it applied first.

pub(crate) use super::error::ReadSolutionsError;

use std::path::Path;

use crate::ndarray::{prelude::*, Array3};
use marlu::{
    fitsio::{errors::Error as FitsioError, hdu::FitsHdu, FitsFile},
    hifitime::Epoch,
    num_complex::Complex,
    Jones,
};

/// The corrections that hyperdrive applied to raw MWA data before it derived
/// the calibration solutions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawDataCorrections {
    /// hyperdrive's name for the polyphase filter bank passband gains it
    /// corrected, e.g. `none`, `jake`, `jake_oversampled`, `cotter2014`
    pub pfb_flavour: String,
    /// Whether the digital gains were corrected
    pub digital_gains: bool,
    /// Whether the cable lengths were corrected
    pub cable_length: bool,
    /// Whether the geometric delays were corrected
    pub geometric: bool,
}

/// The relevant information contained within a hyperdrive fits calibration
/// solutions file.
#[derive(Clone, Debug)]
pub struct HyperdriveSols {
    /// a three dimensional array of jones matrix calibration solutions with
    /// dimensions `[timeblock][tile][chanblock]`
    pub di_jones: Array3<Jones<f64>>,

    /// The observation that the solutions were derived from
    pub obsid: Option<u32>,

    /// The centre frequency \[Hz\] of each chanblock, if the file has them.
    pub chanblock_freqs_hz: Option<Vec<f64>>,

    /// The name of each tile, if the file has them.
    pub tile_names: Option<Vec<String>>,

    /// The start time of each timeblock, if the file has them.
    pub start_timestamps: Vec<Epoch>,

    /// The corrections applied before calibration, if the file has them.
    pub raw_data_corrections: Option<RawDataCorrections>,
}

/// Is this the error cfitsio gives for something that isn't in the file?
fn is_missing(e: &FitsioError) -> bool {
    match e {
        // 202: keyword not found, 219: column not found, 301: HDU not found
        FitsioError::Fits(fe) => matches!(fe.status, 202 | 219 | 301),
        FitsioError::Message(m) => m.starts_with("Cannot find"),
        _ => false,
    }
}

/// Turn "not in the file" into `None`, keep any other error.
fn optional<T>(result: Result<T, FitsioError>) -> Result<Option<T>, FitsioError> {
    match result {
        Ok(v) => Ok(Some(v)),
        Err(e) if is_missing(&e) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Read a column from a table that might not be there.
fn optional_col<T: marlu::fitsio::tables::ReadsCol>(
    fptr: &mut FitsFile,
    hdu_name: &str,
    col_name: &str,
) -> Result<Option<Vec<T>>, FitsioError> {
    let Some(hdu): Option<FitsHdu> = optional(fptr.hdu(hdu_name))? else {
        return Ok(None);
    };
    optional(hdu.read_col::<T>(fptr, col_name))
}

impl HyperdriveSols {
    /// Read a hyperdrive fits calibration solutions file.
    ///
    /// # Errors
    ///
    /// Can throw [`ReadSolutionsError`] if the file can't be read as fits, has
    /// no `SOLUTIONS` image, or its tables don't match the shape of that image.
    pub fn read<T: AsRef<Path>>(file: T) -> Result<Self, ReadSolutionsError> {
        let file_str = file.as_ref().display().to_string();
        let mut fptr = FitsFile::open(&file)?;

        let hdu = fptr.primary_hdu()?;
        let obsid = optional(hdu.read_key::<String>(&mut fptr, "OBSID"))?
            .and_then(|s| s.trim().parse::<u32>().ok());
        let yes = |v: &str| v.trim() == "Y";
        let raw_data_corrections = match (
            optional(hdu.read_key::<String>(&mut fptr, "PFB"))?,
            optional(hdu.read_key::<String>(&mut fptr, "D_GAINS"))?,
            optional(hdu.read_key::<String>(&mut fptr, "CABLELEN"))?,
            optional(hdu.read_key::<String>(&mut fptr, "GEOMETRY"))?,
        ) {
            (Some(pfb), Some(d), Some(c), Some(g)) => Some(RawDataCorrections {
                pfb_flavour: pfb.trim().to_lowercase(),
                digital_gains: yes(&d),
                cable_length: yes(&c),
                geometric: yes(&g),
            }),
            _ => None,
        };

        let hdu = fptr.hdu("SOLUTIONS")?;
        let mut naxis = [0_usize; 4];
        for (i, n) in naxis.iter_mut().enumerate() {
            *n = hdu.read_key::<i64>(&mut fptr, &format!("NAXIS{}", i + 1))? as usize;
        }
        let [num_floats, num_chanblocks, num_tiles, num_timeblocks] = naxis;
        if num_floats != 8 {
            return Err(ReadSolutionsError::BadShape {
                file: file_str,
                thing: "the number of floats in each Jones matrix (NAXIS1 of SOLUTIONS)",
                expected: 8,
                received: num_floats,
            });
        }
        let di_jones_vec: Vec<f64> = hdu.read_image(&mut fptr)?;
        let di_jones = Array4::from_shape_vec(
            (num_timeblocks, num_tiles, num_chanblocks, num_floats),
            di_jones_vec,
        )
        .expect("cfitsio reads the whole image")
        .map_axis(Axis(3), |view| {
            Jones::from([
                Complex::new(view[0], view[1]),
                Complex::new(view[2], view[3]),
                Complex::new(view[4], view[5]),
                Complex::new(view[6], view[7]),
            ])
        });

        // Like hyperdrive, we say there are no frequencies if any are NaN.
        let chanblock_freqs_hz = optional_col::<f64>(&mut fptr, "CHANBLOCKS", "Freq")?
            .filter(|freqs| !freqs.iter().any(|f| f.is_nan()));
        if let Some(freqs) = chanblock_freqs_hz.as_ref() {
            if freqs.len() != num_chanblocks {
                return Err(ReadSolutionsError::BadShape {
                    file: file_str,
                    thing: "the number of chanblock frequencies in CHANBLOCKS",
                    expected: num_chanblocks,
                    received: freqs.len(),
                });
            }
        }

        let tile_names = optional_col::<String>(&mut fptr, "TILES", "TileName")?;
        if let Some(names) = tile_names.as_ref() {
            if names.len() != num_tiles {
                return Err(ReadSolutionsError::BadShape {
                    file: file_str,
                    thing: "the number of tiles in TILES",
                    expected: num_tiles,
                    received: names.len(),
                });
            }
        }

        // all zeros means the times are not known.
        let start_timestamps = optional_col::<f64>(&mut fptr, "TIMEBLOCKS", "Start")?
            .filter(|starts| !starts.iter().all(|&t| t == 0.0))
            .map(|starts| starts.into_iter().map(Epoch::from_gpst_seconds).collect())
            .unwrap_or_default();

        Ok(Self {
            di_jones,
            obsid,
            chanblock_freqs_hz,
            tile_names,
            start_timestamps,
            raw_data_corrections,
        })
    }
}

impl HyperdriveSols {
    /// Pick out the solutions for the given tiles and frequencies, as an array
    /// with dimensions `[tile][channel]` that
    /// [`crate::calibration::apply_di_calsol`] can apply to data with
    /// `data_freqs_hz.len()` channels.
    ///
    /// Tiles are matched by name if the file has tile names, and by index
    /// otherwise. Channels are matched by frequency if the file has chanblock
    /// frequencies (see [`match_chanblocks`]), otherwise the chanblocks are
    /// assumed to span the same band as the data. Tiles and channels that the
    /// file has no solutions for are NaN, which flags them.
    ///
    /// The number of channels in the result divides `data_freqs_hz.len()`. It is
    /// smaller when each chanblock spans the same whole number of data channels.
    ///
    /// # Errors
    ///
    /// Can throw [`ReadSolutionsError`] if the file has more than one timeblock,
    /// or its solutions can't be matched to the data.
    pub fn select(
        &self,
        tile_names: &[String],
        data_freqs_hz: &[f64],
    ) -> Result<Array2<Jones<f64>>, ReadSolutionsError> {
        let (num_timeblocks, num_tiles, num_chanblocks) = self.di_jones.dim();
        if num_timeblocks != 1 {
            return Err(ReadSolutionsError::BadSelection {
                reason: format!(
                    "the solutions have {num_timeblocks} timeblocks, but Birli can only apply \
                    solutions with one timeblock"
                ),
            });
        }
        let sols = self.di_jones.index_axis(Axis(0), 0);

        // for each tile in the data, the index of its solutions
        let tile_idxs: Vec<Option<usize>> = match self.tile_names.as_ref() {
            Some(sol_names) => tile_names
                .iter()
                .map(|name| sol_names.iter().position(|n| n.trim() == name.trim()))
                .collect(),
            None => {
                if num_tiles != tile_names.len() {
                    return Err(ReadSolutionsError::BadSelection {
                        reason: format!(
                            "the solutions have {num_tiles} tiles and no tile names, but the \
                            observation has {} tiles",
                            tile_names.len()
                        ),
                    });
                }
                (0..num_tiles).map(Some).collect()
            }
        };
        let missing_tiles: Vec<&str> = tile_names
            .iter()
            .zip(tile_idxs.iter())
            .filter_map(|(name, idx)| idx.is_none().then_some(name.as_str()))
            .collect();
        if missing_tiles.len() == tile_names.len() {
            return Err(ReadSolutionsError::BadSelection {
                reason: "none of the tiles in the observation are in the solutions".into(),
            });
        }
        if !missing_tiles.is_empty() {
            log::warn!(
                "These tiles are not in the calibration solutions and will be flagged: {}",
                missing_tiles.join(", ")
            );
        }

        // for each channel of the result, the index of its chanblock
        let chan_idxs: Vec<Option<usize>> = match self.chanblock_freqs_hz.as_ref() {
            Some(chanblock_freqs_hz) => {
                let matches = match_chanblocks(chanblock_freqs_hz, data_freqs_hz);
                if matches.iter().all(Option::is_none) {
                    return Err(ReadSolutionsError::BadSelection {
                        reason: format!(
                            "the solutions cover {:.3} to {:.3} MHz, which does not overlap \
                            with the selected data, {:.3} to {:.3} MHz",
                            chanblock_freqs_hz.first().copied().unwrap_or(f64::NAN) / 1e6,
                            chanblock_freqs_hz.last().copied().unwrap_or(f64::NAN) / 1e6,
                            data_freqs_hz.first().copied().unwrap_or(f64::NAN) / 1e6,
                            data_freqs_hz.last().copied().unwrap_or(f64::NAN) / 1e6,
                        ),
                    });
                }
                let num_missing = matches.iter().filter(|m| m.is_none()).count();
                if num_missing > 0 {
                    log::warn!(
                        "{num_missing} of {} selected channels are outside the frequencies of \
                        the calibration solutions and will be flagged",
                        matches.len()
                    );
                }
                let ratio = regular_chanblock_ratio(&matches).unwrap_or(1);
                matches.into_iter().step_by(ratio).collect()
            }
            None => {
                if num_chanblocks == 0 || !data_freqs_hz.len().is_multiple_of(num_chanblocks) {
                    return Err(ReadSolutionsError::BadSelection {
                        reason: format!(
                            "the solutions have {num_chanblocks} chanblocks and no chanblock \
                            frequencies, which does not divide the {} selected channels",
                            data_freqs_hz.len()
                        ),
                    });
                }
                (0..num_chanblocks).map(Some).collect()
            }
        };

        let nan = Jones::from([Complex::new(f64::NAN, f64::NAN); 4]);
        Ok(Array2::from_shape_fn(
            (tile_idxs.len(), chan_idxs.len()),
            |(tile, chan)| match (tile_idxs[tile], chan_idxs[chan]) {
                (Some(tile_idx), Some(chan_idx)) => sols[(tile_idx, chan_idx)],
                _ => nan,
            },
        ))
    }
}

/// For each of `data_freqs_hz`, the index of the chanblock that contains it, or
/// `None` if no chanblock does.
///
/// A chanblock contains the frequencies within half of the chanblock spacing of
/// its centre. If there is only one chanblock its width isn't known, so it
/// contains everything.
///
/// `chanblock_freqs_hz` must be sorted in ascending order.
pub fn match_chanblocks(chanblock_freqs_hz: &[f64], data_freqs_hz: &[f64]) -> Vec<Option<usize>> {
    let spacing = chanblock_freqs_hz
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f64::INFINITY, f64::min);
    // allow for a little rounding in the frequencies
    let tolerance = spacing * (0.5 + 1e-6);
    data_freqs_hz
        .iter()
        .map(|&freq| {
            let right = chanblock_freqs_hz.partition_point(|&f| f < freq);
            let nearest = [right.checked_sub(1), Some(right)]
                .into_iter()
                .flatten()
                .filter(|&i| i < chanblock_freqs_hz.len())
                .min_by(|&a, &b| {
                    (chanblock_freqs_hz[a] - freq)
                        .abs()
                        .total_cmp(&(chanblock_freqs_hz[b] - freq).abs())
                })?;
            ((chanblock_freqs_hz[nearest] - freq).abs() <= tolerance).then_some(nearest)
        })
        .collect()
}

/// If the channels in `matches` (see [`match_chanblocks`]) come in runs of
/// equal length, each with a single chanblock, return that length. This is when
/// the solutions can be applied without expanding them to the resolution of the
/// data.
pub fn regular_chanblock_ratio(matches: &[Option<usize>]) -> Option<usize> {
    // the length of the first run of channels that share a chanblock
    let start = matches.iter().position(Option::is_some)?;
    let run = matches[start..]
        .iter()
        .position(|m| *m != matches[start])
        .unwrap_or(matches.len() - start);
    if !matches.len().is_multiple_of(run) {
        return None;
    }
    let mut previous = None;
    for chunk in matches.chunks(run) {
        if chunk.iter().any(|m| *m != chunk[0]) {
            return None;
        }
        // neighbouring runs must not share a chanblock, unless both are
        // outside the solutions.
        if chunk[0].is_some() && previous == Some(chunk[0]) {
            return None;
        }
        previous = Some(chunk[0]);
    }
    Some(run)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use marlu::fitsio::{
        images::{ImageDescription, ImageType},
        tables::{ColumnDataType, ColumnDescription},
    };
    use tempfile::tempdir;

    /// Write a hyperdrive-style fits solutions file, for testing.
    pub(crate) fn write_hyperdrive_fits(
        path: &Path,
        di_jones: ArrayView3<Jones<f64>>,
        chanblock_freqs_hz: Option<&[f64]>,
        tile_names: Option<&[String]>,
        raw_data_corrections: Option<&RawDataCorrections>,
        obsid: Option<u32>,
    ) {
        let mut fptr = FitsFile::create(path).overwrite().open().unwrap();
        let hdu = fptr.primary_hdu().unwrap();
        if let Some(obsid) = obsid {
            hdu.write_key(&mut fptr, "OBSID", format!("{obsid}"))
                .unwrap();
        }
        if let Some(c) = raw_data_corrections {
            let yn = |b: bool| if b { "Y" } else { "N" };
            hdu.write_key(&mut fptr, "PFB", c.pfb_flavour.as_str())
                .unwrap();
            hdu.write_key(&mut fptr, "D_GAINS", yn(c.digital_gains))
                .unwrap();
            hdu.write_key(&mut fptr, "CABLELEN", yn(c.cable_length))
                .unwrap();
            hdu.write_key(&mut fptr, "GEOMETRY", yn(c.geometric))
                .unwrap();
        }
        let (num_timeblocks, num_tiles, num_chanblocks) = di_jones.dim();
        let dim = [num_timeblocks, num_tiles, num_chanblocks, 8];
        let image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &dim,
        };
        let hdu = fptr.create_image("SOLUTIONS", &image_description).unwrap();
        let floats: Vec<f64> = di_jones
            .iter()
            .flat_map(|j| j.iter().flat_map(|c| [c.re, c.im]).collect::<Vec<_>>())
            .collect();
        hdu.write_image(&mut fptr, &floats).unwrap();

        if let Some(names) = tile_names {
            let cols = [
                ColumnDescription::new("Antenna")
                    .with_type(ColumnDataType::Int)
                    .create()
                    .unwrap(),
                ColumnDescription::new("TileName")
                    .with_type(ColumnDataType::String)
                    .that_repeats(8)
                    .create()
                    .unwrap(),
            ];
            let hdu = fptr.create_table("TILES", &cols).unwrap();
            let idxs: Vec<i32> = (0..names.len() as i32).collect();
            hdu.write_col(&mut fptr, "Antenna", &idxs).unwrap();
            hdu.write_col(&mut fptr, "TileName", names).unwrap();
        }
        if let Some(freqs) = chanblock_freqs_hz {
            let cols = [
                ColumnDescription::new("Index")
                    .with_type(ColumnDataType::Int)
                    .create()
                    .unwrap(),
                ColumnDescription::new("Freq")
                    .with_type(ColumnDataType::Double)
                    .create()
                    .unwrap(),
            ];
            let hdu = fptr.create_table("CHANBLOCKS", &cols).unwrap();
            let idxs: Vec<i32> = (0..freqs.len() as i32).collect();
            hdu.write_col(&mut fptr, "Index", &idxs).unwrap();
            hdu.write_col(&mut fptr, "Freq", freqs).unwrap();
        }
    }

    fn test_jones(num_timeblocks: usize, num_tiles: usize, num_chans: usize) -> Array3<Jones<f64>> {
        Array3::from_shape_fn((num_timeblocks, num_tiles, num_chans), |(t, a, c)| {
            let v = (t * 100 + a * 10 + c) as f64;
            Jones::from([
                Complex::new(v, 1.0),
                Complex::new(v, 2.0),
                Complex::new(v, 3.0),
                Complex::new(v, 4.0),
            ])
        })
    }

    #[test]
    fn test_read_hyperdrive_fits_full() {
        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("sols.fits");
        let mut di_jones = test_jones(1, 3, 4);
        // a flagged tile
        di_jones
            .slice_mut(s![.., 1, ..])
            .fill(Jones::from([Complex::new(f64::NAN, f64::NAN); 4]));
        let freqs = [150e6, 150.04e6, 150.08e6, 150.12e6];
        let names: Vec<String> = ["Tile011", "Tile012", "Tile013"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let corrections = RawDataCorrections {
            pfb_flavour: "jake".into(),
            digital_gains: true,
            cable_length: false,
            geometric: true,
        };
        write_hyperdrive_fits(
            &path,
            di_jones.view(),
            Some(&freqs),
            Some(&names),
            Some(&corrections),
            Some(1_254_670_392),
        );

        let sols = HyperdriveSols::read(&path).unwrap();
        assert_eq!(sols.di_jones.dim(), (1, 3, 4));
        assert_eq!(sols.di_jones[(0, 2, 3)], di_jones[(0, 2, 3)]);
        assert_eq!(sols.di_jones[(0, 0, 1)][2], Complex::new(1.0, 3.0));
        assert!(sols.di_jones[(0, 1, 0)].any_nan());
        assert_eq!(sols.obsid, Some(1_254_670_392));
        assert_eq!(sols.chanblock_freqs_hz.as_deref(), Some(&freqs[..]));
        assert_eq!(sols.tile_names, Some(names));
        assert_eq!(sols.raw_data_corrections, Some(corrections));
        assert!(sols.start_timestamps.is_empty());
    }

    /// Timeblock start times are read as GPS seconds, and a table of zeros
    /// means the times are not known.
    #[test]
    fn test_read_hyperdrive_fits_timeblocks() {
        let tmp_dir = tempdir().unwrap();
        let di_jones = test_jones(2, 2, 3);
        for (name, starts, expected) in [
            (
                "known.fits",
                [1_254_670_392.0, 1_254_670_400.0],
                vec![
                    Epoch::from_gpst_seconds(1_254_670_392.0),
                    Epoch::from_gpst_seconds(1_254_670_400.0),
                ],
            ),
            ("unknown.fits", [0.0, 0.0], vec![]),
        ] {
            let path = tmp_dir.path().join(name);
            write_hyperdrive_fits(&path, di_jones.view(), None, None, None, None);
            let mut fptr = FitsFile::edit(&path).unwrap();
            let cols = [ColumnDescription::new("Start")
                .with_type(ColumnDataType::Double)
                .create()
                .unwrap()];
            let hdu = fptr.create_table("TIMEBLOCKS", &cols).unwrap();
            hdu.write_col(&mut fptr, "Start", &starts).unwrap();
            drop(fptr);

            let sols = HyperdriveSols::read(&path).unwrap();
            assert_eq!(sols.start_timestamps, expected, "{name}");
        }
    }

    #[test]
    fn test_read_hyperdrive_fits_minimal() {
        let tmp_dir = tempdir().unwrap();
        let path = tmp_dir.path().join("sols.fits");
        let di_jones = test_jones(2, 2, 3);
        write_hyperdrive_fits(&path, di_jones.view(), None, None, None, None);

        let sols = HyperdriveSols::read(&path).unwrap();
        assert_eq!(sols.di_jones, di_jones);
        assert_eq!(sols.obsid, None);
        assert_eq!(sols.chanblock_freqs_hz, None);
        assert_eq!(sols.tile_names, None);
        assert_eq!(sols.raw_data_corrections, None);
    }

    #[test]
    fn test_read_hyperdrive_fits_bad_files() {
        let tmp_dir = tempdir().unwrap();
        // not a fits file
        let path = tmp_dir.path().join("not.fits");
        std::fs::write(&path, b"not a fits file").unwrap();
        assert!(matches!(
            HyperdriveSols::read(&path),
            Err(ReadSolutionsError::Fitsio(_))
        ));
        // more frequencies than chanblocks
        let path = tmp_dir.path().join("sols.fits");
        let di_jones = test_jones(1, 2, 3);
        write_hyperdrive_fits(&path, di_jones.view(), Some(&[1e6, 2e6]), None, None, None);
        assert!(matches!(
            HyperdriveSols::read(&path),
            Err(ReadSolutionsError::BadShape { .. })
        ));
    }

    #[test]
    fn test_select_by_name_and_frequency() {
        let di_jones = test_jones(1, 3, 3);
        let names = |v: &[&str]| v.iter().map(|s| (*s).to_string()).collect::<Vec<_>>();
        let sols = HyperdriveSols {
            di_jones: di_jones.clone(),
            obsid: None,
            chanblock_freqs_hz: Some(vec![100.0e3, 140.0e3, 180.0e3]),
            // fits pads strings with spaces
            tile_names: Some(names(&["Tile011 ", "Tile012", "Tile013"])),
            start_timestamps: vec![],
            raw_data_corrections: None,
        };

        // the data has the tiles in a different order, a tile that is not in the
        // solutions, and 10 kHz channels from the middle of the first chanblock
        // to past the last one.
        let data_names = names(&["Tile013", "Tile099", "Tile011"]);
        let data_freqs: Vec<f64> = (0..12).map(|i| 105.0e3 + 10.0e3 * i as f64).collect();
        let selected = sols.select(&data_names, &data_freqs).unwrap();
        // not regular, so one solution for each data channel
        assert_eq!(selected.dim(), (3, 12));
        assert_eq!(selected[(0, 0)], di_jones[(0, 2, 0)]);
        assert_eq!(selected[(0, 1)], di_jones[(0, 2, 0)]);
        assert_eq!(selected[(0, 2)], di_jones[(0, 2, 1)]);
        assert_eq!(selected[(2, 9)], di_jones[(0, 0, 2)]);
        assert!(selected[(2, 10)].any_nan());
        assert!(selected.row(1).iter().all(|j| j.any_nan()));

        // data aligned with the chanblocks: one solution for each chanblock
        let data_freqs: Vec<f64> = (0..8).map(|i| 125.0e3 + 10.0e3 * i as f64).collect();
        let selected = sols.select(&data_names, &data_freqs).unwrap();
        assert_eq!(selected.dim(), (3, 2));
        assert_eq!(selected[(2, 0)], di_jones[(0, 0, 1)]);
        assert_eq!(selected[(2, 1)], di_jones[(0, 0, 2)]);

        // no overlap in frequency, or in tiles
        assert!(matches!(
            sols.select(&data_names, &[1e9, 2e9]),
            Err(ReadSolutionsError::BadSelection { .. })
        ));
        assert!(matches!(
            sols.select(&names(&["Tile098", "Tile099"]), &data_freqs),
            Err(ReadSolutionsError::BadSelection { .. })
        ));
    }

    #[test]
    fn test_select_without_names_or_frequencies() {
        let names = |n: usize| (0..n).map(|i| format!("Tile{i:03}")).collect::<Vec<_>>();
        let mut sols = HyperdriveSols {
            di_jones: test_jones(1, 3, 2),
            obsid: None,
            chanblock_freqs_hz: None,
            tile_names: None,
            start_timestamps: vec![],
            raw_data_corrections: None,
        };
        let selected = sols.select(&names(3), &[1e6, 2e6, 3e6, 4e6]).unwrap();
        assert_eq!(selected, sols.di_jones.index_axis(Axis(0), 0));
        // wrong number of tiles
        assert!(sols.select(&names(4), &[1e6, 2e6]).is_err());
        // 2 chanblocks does not divide 3 channels
        assert!(sols.select(&names(3), &[1e6, 2e6, 3e6]).is_err());
        // more than one timeblock
        sols.di_jones = test_jones(2, 3, 2);
        assert!(sols.select(&names(3), &[1e6, 2e6]).is_err());
    }

    #[test]
    fn test_missing_or_broken() {
        use marlu::fitsio::errors::FitsError;
        let fits = |status| {
            FitsioError::Fits(FitsError {
                status,
                message: String::new(),
            })
        };
        assert!(is_missing(&fits(202)));
        assert!(is_missing(&fits(301)));
        assert!(!is_missing(&fits(104)));
        assert!(is_missing(&FitsioError::Message(
            "Cannot find column \"Freq\"".into()
        )));
        assert!(!is_missing(&FitsioError::Message("something else".into())));
        assert!(!is_missing(&FitsioError::UnlockError));

        assert_eq!(optional(Ok(1)).unwrap(), Some(1));
        assert_eq!(optional::<i32>(Err(fits(202))).unwrap(), None);
        assert!(optional::<i32>(Err(fits(104))).is_err());
    }

    #[test]
    fn test_read_hyperdrive_fits_bad_shapes() {
        let tmp_dir = tempdir().unwrap();
        // more tile names than tiles
        let path = tmp_dir.path().join("tiles.fits");
        let names: Vec<String> = (0..3).map(|i| format!("Tile{i:03}")).collect();
        write_hyperdrive_fits(
            &path,
            test_jones(1, 2, 3).view(),
            None,
            Some(&names),
            None,
            None,
        );
        assert!(matches!(
            HyperdriveSols::read(&path),
            Err(ReadSolutionsError::BadShape { .. })
        ));

        // not 8 floats for each Jones matrix
        let path = tmp_dir.path().join("floats.fits");
        let mut fptr = FitsFile::create(&path).overwrite().open().unwrap();
        let image_description = ImageDescription {
            data_type: ImageType::Double,
            dimensions: &[1, 2, 3, 4],
        };
        let hdu = fptr.create_image("SOLUTIONS", &image_description).unwrap();
        hdu.write_image(&mut fptr, &[0.0_f64; 24]).unwrap();
        drop(fptr);
        let err = HyperdriveSols::read(&path).unwrap_err();
        assert!(format!("{err}").contains("NAXIS1 of SOLUTIONS"), "{err}");
    }

    #[test]
    fn test_match_chanblocks() {
        // 40 kHz chanblocks, 10 kHz data that starts before and ends after them
        let chanblocks = [100.0e3, 140.0e3, 180.0e3];
        let data: Vec<f64> = (0..16).map(|i| 75.0e3 + 10.0e3 * i as f64).collect();
        let matches = match_chanblocks(&chanblocks, &data);
        #[rustfmt::skip]
        assert_eq!(matches, vec![
            None,
            Some(0), Some(0), Some(0), Some(0),
            Some(1), Some(1), Some(1), Some(1),
            Some(2), Some(2), Some(2), Some(2),
            None, None, None,
        ]);
        assert_eq!(regular_chanblock_ratio(&matches), None);
        assert_eq!(regular_chanblock_ratio(&matches[1..13]), Some(4));
        // same resolution
        let matches = match_chanblocks(&chanblocks, &chanblocks);
        assert_eq!(matches, vec![Some(0), Some(1), Some(2)]);
        assert_eq!(regular_chanblock_ratio(&matches), Some(1));
        // a single chanblock has no known width
        assert_eq!(match_chanblocks(&[1e6], &[0.0, 2e9]), vec![Some(0); 2]);
        assert_eq!(regular_chanblock_ratio(&[Some(0); 2]), Some(2));
        // two runs with the same chanblock can't be told apart
        assert_eq!(
            regular_chanblock_ratio(&[Some(0), Some(1), Some(1), Some(2)]),
            None
        );
        // a run that the channels don't divide into
        assert_eq!(regular_chanblock_ratio(&[Some(0), Some(0), Some(1)]), None);
        // runs outside the solutions are fine
        assert_eq!(
            regular_chanblock_ratio(&[None, None, None, None, Some(0), Some(0)]),
            Some(2)
        );
    }
}
