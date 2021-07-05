// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Module for uvfits file format reading and writing
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

use std::collections::HashSet;
use std::ffi::CString;
use std::path::Path;

use erfa_sys::{ERFA_DJM0, ERFA_WGS84};
use fitsio::{errors::check_status as fits_check_status, FitsFile};
use hifitime::Epoch;
use ndarray::prelude::*;

use super::error::UvfitsWriteError;
use crate::coord::ENH;
use crate::{
    coord::{RADec, XyzGeocentric, XyzGeodetic, UVW},
    math::{cross_correlation_baseline_to_tiles, num_tiles_from_num_baselines},
};

// use mwa_hyperdrive_core::{
//     Jones, RADec, XyzGeocentric, XyzGeodetic, UVW,
// };
use mwalib::{fitsio, fitsio_sys, Antenna, CorrelatorContext, SPEED_OF_LIGHT_IN_VACUUM_M_PER_S};

/// From a `hifitime` [Epoch], get a formatted date string with the hours,
/// minutes and seconds set to 0.
fn get_truncated_date_string(epoch: &Epoch) -> String {
    let (year, month, day, _, _, _, _) = epoch.as_gregorian_utc();
    format!(
        "{year}-{month}-{day}T00:00:00.0",
        year = year,
        month = month,
        day = day
    )
}

/// Helper function to convert strings into pointers of C strings.
fn rust_strings_to_c_strings<T: AsRef<str>>(
    strings: &[T],
) -> Result<Vec<*mut i8>, std::ffi::NulError> {
    let mut c_strings = Vec::with_capacity(strings.len());
    for s in strings {
        let rust_str = s.as_ref().to_owned();
        let c_str = CString::new(rust_str)?;
        c_strings.push(c_str.into_raw());
    }
    Ok(c_strings)
}

/// Encode a baseline into the uvfits format. Use the miriad convention to
/// handle more than 255 antennas (up to 2048). This is backwards compatible
/// with the standard UVFITS convention. Antenna indices start at 1.
// Shamelessly copied from the RTS, originally written by Randall Wayth.
fn encode_uvfits_baseline(ant1: usize, ant2: usize) -> usize {
    if ant2 > 255 {
        ant1 * 2048 + ant2 + 65_536
    } else {
        ant1 * 256 + ant2
    }
}

/// Decode a uvfits baseline into the antennas that formed it. Antenna indices
/// start at 1.
fn decode_uvfits_baseline(bl: usize) -> (usize, usize) {
    if bl < 65_535 {
        let ant2 = bl % 256;
        let ant1 = (bl - ant2) / 256;
        (ant1, ant2)
    } else {
        let ant2 = (bl - 65_536) % 2048;
        let ant1 = (bl - ant2 - 65_536) / 2048;
        (ant1, ant2)
    }
}

/// A helper struct to write out a uvfits file.
///
/// TODO: make writer and reader a single class?
pub(crate) struct UvfitsWriter<'a> {
    /// The path to the uvifts file.
    path: &'a Path,

    /// The number of timesteps to be written.    
    num_timesteps: usize,

    /// The number of baselines per timestep to be written.
    num_baselines: usize,

    /// The number of fine channels per baseline to be written. uvfits has no
    /// notion of "fine channel" or "coarse channel".
    num_chans: usize,

    /// The number of uvfits rows. This is equal to `num_timesteps` *
    /// `num_baselines`.
    total_num_rows: usize,

    /// The number of uvfits rows that have currently been written.
    current_num_rows: usize,

    /// The frequency at the centre of the bandwidth \[Hz\].
    centre_freq: f64,

    /// A `hifitime` [Epoch] struct associated with the first timestep of the
    /// data.
    start_epoch: &'a Epoch,
}

impl<'a> UvfitsWriter<'a> {
    /// Create a new uvfits file at the specified filename.
    pub(crate) fn new(
        filename: &'a Path,
        num_timesteps: usize,
        num_baselines: usize,
        num_chans: usize,
        start_epoch: &'a Epoch,
        fine_chan_width_hz: f64,
        centre_freq_hz: f64,
        centre_freq_chan: usize,
        phase_centre: &RADec,
        obs_name: Option<&str>,
    ) -> Result<Self, UvfitsWriteError> {
        // Delete any file that already exists.
        if filename.exists() {
            std::fs::remove_file(&filename)?;
        }

        // Create a new fits file.
        let mut status = 0;
        let c_filename = CString::new(filename.to_str().unwrap())?;
        let mut fptr = std::ptr::null_mut();
        unsafe {
            fitsio_sys::ffinit(
                &mut fptr as *mut *mut _, /* O - FITS file pointer                   */
                c_filename.as_ptr(),      /* I - name of file to create              */
                &mut status,              /* IO - error status                       */
            );
        }
        fits_check_status(status)?;

        // Initialise the group header. Copied from cotter. -32 means FLOAT_IMG.
        let naxis = 6;
        let mut naxes = [0, 3, 4, num_chans as i64, 1, 1];
        let num_group_params = 5;
        let total_num_rows = num_timesteps * num_baselines;
        unsafe {
            fitsio_sys::ffphpr(
                fptr,                  /* I - FITS file pointer                        */
                1,                     /* I - does file conform to FITS standard? 1/0  */
                -32,                   /* I - number of bits per data value pixel      */
                naxis,                 /* I - number of axes in the data array         */
                naxes.as_mut_ptr(),    /* I - length of each data axis                 */
                num_group_params,      /* I - number of group parameters (usually 0)   */
                total_num_rows as i64, /* I - number of random groups (usually 1 or 0) */
                1,                     /* I - may FITS file have extensions?           */
                &mut status,           /* IO - error status                            */
            );
        }
        fits_check_status(status)?;

        // Finally close the file.
        unsafe {
            fitsio_sys::ffclos(fptr, &mut status);
        }
        fits_check_status(status)?;

        // Open the fits file with rust-fitsio.
        let mut u = FitsFile::edit(&filename)?;
        let hdu = u.hdu(0)?;
        hdu.write_key(&mut u, "BSCALE", 1.0)?;

        // Set header names and scales.
        for (i, &param) in ["UU", "VV", "WW", "BASELINE", "DATE"].iter().enumerate() {
            let ii = i + 1;
            hdu.write_key(&mut u, &format!("PTYPE{}", ii), param)?;
            hdu.write_key(&mut u, &format!("PSCAL{}", ii), 1.0)?;
            if param != "DATE" {
                hdu.write_key(&mut u, &format!("PZERO{}", ii), 0.0)?;
            } else {
                // Set the zero level for the DATE column.
                hdu.write_key(
                    &mut u,
                    &format!("PZERO{}", ii),
                    start_epoch.as_jde_utc_days().floor() + 0.5,
                )?;
            }
        }
        hdu.write_key(&mut u, "DATE-OBS", get_truncated_date_string(&start_epoch))?;

        // Dimensions.
        hdu.write_key(&mut u, "CTYPE2", "COMPLEX")?;
        hdu.write_key(&mut u, "CRVAL2", 1.0)?;
        hdu.write_key(&mut u, "CRPIX2", 1.0)?;
        hdu.write_key(&mut u, "CDELT2", 1.0)?;

        // Linearly polarised.
        hdu.write_key(&mut u, "CTYPE3", "STOKES")?;
        hdu.write_key(&mut u, "CRVAL3", -5)?;
        hdu.write_key(&mut u, "CDELT3", -1)?;
        hdu.write_key(&mut u, "CRPIX3", 1.0)?;

        hdu.write_key(&mut u, "CTYPE4", "FREQ")?;
        hdu.write_key(&mut u, "CRVAL4", centre_freq_hz)?;
        hdu.write_key(&mut u, "CDELT4", fine_chan_width_hz)?;
        hdu.write_key(&mut u, "CRPIX4", centre_freq_chan as u64 + 1)?;

        hdu.write_key(&mut u, "CTYPE5", "RA")?;
        hdu.write_key(&mut u, "CRVAL5", phase_centre.ra.to_degrees())?;
        hdu.write_key(&mut u, "CDELT5", 1)?;
        hdu.write_key(&mut u, "CRPIX5", 1)?;

        hdu.write_key(&mut u, "CTYPE6", "DEC")?;
        hdu.write_key(&mut u, "CRVAL6", phase_centre.dec.to_degrees())?;
        hdu.write_key(&mut u, "CDELT6", 1)?;
        hdu.write_key(&mut u, "CRPIX6", 1)?;

        hdu.write_key(&mut u, "OBSRA", phase_centre.ra.to_degrees())?;
        hdu.write_key(&mut u, "OBSDEC", phase_centre.dec.to_degrees())?;
        hdu.write_key(&mut u, "EPOCH", 2000.0)?;

        hdu.write_key(&mut u, "OBJECT", obs_name.unwrap_or("Undefined"))?;
        hdu.write_key(&mut u, "TELESCOP", "MWA")?;
        hdu.write_key(&mut u, "INSTRUME", "MWA")?;

        // This is apparently required...
        let history = CString::new("AIPS WTSCAL =  1.0").unwrap();
        unsafe {
            fitsio_sys::ffphis(
                u.as_raw(),       /* I - FITS file pointer  */
                history.as_ptr(), /* I - history string     */
                &mut status,      /* IO - error status      */
            );
        }
        fits_check_status(status)?;

        // Add in version information
        let comment = CString::new(format!(
            "Created by {} v{}",
            env!("CARGO_PKG_NAME"),
            env!("CARGO_PKG_VERSION")
        ))
        .unwrap();
        unsafe {
            fitsio_sys::ffpcom(
                u.as_raw(),       /* I - FITS file pointer   */
                comment.as_ptr(), /* I - comment string      */
                &mut status,      /* IO - error status       */
            );
        }
        fits_check_status(status)?;

        hdu.write_key(&mut u, "SOFTWARE", env!("CARGO_PKG_NAME"))?;
        hdu.write_key(
            &mut u,
            "GITLABEL",
            format!("v{}", env!("CARGO_PKG_VERSION")),
        )?;

        Ok(Self {
            path: filename,
            num_timesteps,
            num_baselines,
            num_chans,
            total_num_rows,
            current_num_rows: 0,
            centre_freq: centre_freq_hz,
            start_epoch,
        })
    }


    /// Opens the associated uvfits file in edit mode, returning the [FitsFile]
    /// struct.
    pub(crate) fn open(&self) -> Result<FitsFile, fitsio::errors::Error> {
        let mut f = FitsFile::edit(&self.path)?;
        // Ensure HDU 0 is opened.
        f.hdu(0)?;
        Ok(f)
    }

    /// Write the antenna table to a uvfits file. Assumes that the array
    /// location is MWA.
    ///
    /// `centre_freq` is the centre frequency of the coarse band that this
    /// uvfits file pertains to. `positions` are the [XyzGeodetic] coordinates
    /// of the MWA tiles. These positions need to have the MWA's "centre" XYZ
    /// coordinates subtracted to make them local XYZ.
    ///
    /// `Self` must have only have a single HDU when this function is called
    /// (true when using methods only provided by `Self`).
    // Derived from cotter.
    pub(crate) fn write_uvfits_antenna_table<T: AsRef<str>>(
        self,
        antenna_names: &[T],
        positions: &[XyzGeodetic],
    ) -> Result<(), UvfitsWriteError> {
        if self.current_num_rows != self.total_num_rows {
            return Err(UvfitsWriteError::NotEnoughRowsWritten {
                current: self.current_num_rows,
                total: self.total_num_rows,
            });
        }

        let mut uvfits = self.open()?;

        // Stuff that a uvfits file always expects?
        let col_names = [
            "ANNAME", "STABXYZ", "NOSTA", "MNTSTA", "STAXOF", "POLTYA", "POLAA", "POLCALA",
            "POLTYB", "POLAB", "POLCALB",
        ];
        let col_formats = [
            "8A", "3D", "1J", "1J", "1E", "1A", "1E", "3E", "1A", "1E", "3E",
        ];
        let col_units = [
            "", "METERS", "", "", "METERS", "", "DEGREES", "", "", "DEGREES", "",
        ];
        let mut c_col_names = rust_strings_to_c_strings(&col_names)?;
        let mut c_col_formats = rust_strings_to_c_strings(&col_formats)?;
        let mut c_col_units = rust_strings_to_c_strings(&col_units)?;
        let extname = CString::new("AIPS AN").unwrap();

        // ffcrtb creates a new binary table in a new HDU. This should be the second
        // HDU, so there should only be one HDU before this function is called.
        let mut status = 0;
        unsafe {
            // BINARY_TBL is 2.
            fitsio_sys::ffcrtb(
                uvfits.as_raw(),            /* I - FITS file pointer                        */
                2,                          /* I - type of table to create                  */
                0,                          /* I - number of rows in the table              */
                11,                         /* I - number of columns in the table           */
                c_col_names.as_mut_ptr(),   /* I - name of each column                      */
                c_col_formats.as_mut_ptr(), /* I - value of TFORMn keyword for each column  */
                c_col_units.as_mut_ptr(),   /* I - value of TUNITn keyword for each column  */
                extname.as_ptr(),           /* I - value of EXTNAME keyword, if any         */
                &mut status,                /* IO - error status                            */
            );
        }
        fits_check_status(status)?;

        // Open the newly-created HDU.
        let hdu = uvfits.hdu(1)?;

        // Set ARRAYX, Y and Z to the MWA's coordinates in XYZ (geocentric). The
        // results here are slightly different to those given by cotter. This is
        // at least partly due to different constants (the altitude is
        // definitely slightly different), but possibly also because ERFA is
        // more accurate than cotter's "homebrewed" Geodetic2XYZ.
        let mut mwa_xyz: [f64; 3] = [0.0; 3];
        unsafe {
            status = erfa_sys::eraGd2gc(
                ERFA_WGS84 as i32,             // ellipsoid identifier (Note 1)
                mwalib::MWA_LONGITUDE_RADIANS, // longitude (radians, east +ve)
                mwalib::MWA_LATITUDE_RADIANS,  // latitude (geodetic, radians, Note 3)
                mwalib::MWA_ALTITUDE_METRES,   // height above ellipsoid (geodetic, Notes 2,3)
                mwa_xyz.as_mut_ptr(),          // geocentric vector (Note 2)
            );
        }
        if status != 0 {
            return Err(UvfitsWriteError::Erfa {
                source_file: file!(),
                source_line: line!(),
                status,
                function: "eraGd2gc",
            });
        }
        let mwa_xyz = XyzGeocentric {
            x: mwa_xyz[0],
            y: mwa_xyz[1],
            z: mwa_xyz[2],
        };

        hdu.write_key(&mut uvfits, "ARRAYX", mwa_xyz.x)?;
        hdu.write_key(&mut uvfits, "ARRAYY", mwa_xyz.y)?;
        hdu.write_key(&mut uvfits, "ARRAYZ", mwa_xyz.z)?;

        hdu.write_key(&mut uvfits, "FREQ", self.centre_freq)?;

        // Get the Greenwich apparent sidereal time from ERFA.
        let mjd = self.start_epoch.as_mjd_utc_days();
        let gst = unsafe { erfa_sys::eraGst06a(ERFA_DJM0, mjd.floor(), ERFA_DJM0, mjd.floor()) }
            .to_degrees();
        hdu.write_key(&mut uvfits, "GSTIA0", gst)?;
        hdu.write_key(&mut uvfits, "DEGPDY", 3.60985e2)?; // Earth's rotation rate

        let date_truncated = get_truncated_date_string(self.start_epoch);
        hdu.write_key(&mut uvfits, "RDATE", date_truncated)?;

        hdu.write_key(&mut uvfits, "POLARX", 0.0)?;
        hdu.write_key(&mut uvfits, "POLARY", 0.0)?;
        hdu.write_key(&mut uvfits, "UT1UTC", 0.0)?;
        hdu.write_key(&mut uvfits, "DATUTC", 0.0)?;

        hdu.write_key(&mut uvfits, "TIMSYS", "UTC")?;
        hdu.write_key(&mut uvfits, "ARRNAM", "MWA")?;
        hdu.write_key(&mut uvfits, "NUMORB", 0)?; // number of orbital parameters in table
        hdu.write_key(&mut uvfits, "NOPCAL", 3)?; // Nr pol calibration values / IF(N_pcal)
        hdu.write_key(&mut uvfits, "FREQID", -1)?; // Frequency setup number
        hdu.write_key(&mut uvfits, "IATUTC", 33.0)?;

        // Assume the station coordinates are "right handed".
        hdu.write_key(&mut uvfits, "XYZHAND", "RIGHT")?;

        let c_antenna_names = rust_strings_to_c_strings(antenna_names)?;

        // Write to the table row by row.
        for (i, pos) in positions.iter().enumerate() {
            let row = i as i64 + 1;
            unsafe {
                // ANNAME. ffpcls = fits_write_col_str
                fitsio_sys::ffpcls(
                    uvfits.as_raw(),                   /* I - FITS file pointer                       */
                    1,   /* I - number of column to write (1 = 1st col) */
                    row, /* I - first row to write (1 = 1st row)        */
                    1,   /* I - first vector element to write (1 = 1st) */
                    1,   /* I - number of strings to write              */
                    [c_antenna_names[i]].as_mut_ptr(), /* I - array of pointers to strings            */
                    &mut status, /* IO - error status                           */
                );
                fits_check_status(status)?;

                let mut c_xyz = [pos.x, pos.y, pos.z];
                // STABXYZ. ffpcld = fits_write_col_dbl
                fitsio_sys::ffpcld(
                    uvfits.as_raw(),    /* I - FITS file pointer                       */
                    2,                  /* I - number of column to write (1 = 1st col) */
                    row,                /* I - first row to write (1 = 1st row)        */
                    1,                  /* I - first vector element to write (1 = 1st) */
                    3,                  /* I - number of values to write               */
                    c_xyz.as_mut_ptr(), /* I - array of values to write                */
                    &mut status,        /* IO - error status                           */
                );
                fits_check_status(status)?;

                // NOSTA. ffpclk = fits_write_col_int
                fitsio_sys::ffpclk(
                    uvfits.as_raw(),           /* I - FITS file pointer                       */
                    3,                         /* I - number of column to write (1 = 1st col) */
                    row,                       /* I - first row to write (1 = 1st row)        */
                    1,                         /* I - first vector element to write (1 = 1st) */
                    1,                         /* I - number of values to write               */
                    [row as i32].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,               /* IO - error status                           */
                );
                fits_check_status(status)?;

                // MNTSTA
                fitsio_sys::ffpclk(
                    uvfits.as_raw(),  /* I - FITS file pointer                       */
                    4,                /* I - number of column to write (1 = 1st col) */
                    row,              /* I - first row to write (1 = 1st row)        */
                    1,                /* I - first vector element to write (1 = 1st) */
                    1,                /* I - number of values to write               */
                    [0].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,      /* IO - error status                           */
                );
                fits_check_status(status)?;

                // No row 5?
                // POLTYA
                fitsio_sys::ffpcls(
                    uvfits.as_raw(), /* I - FITS file pointer                       */
                    6,               /* I - number of column to write (1 = 1st col) */
                    row,             /* I - first row to write (1 = 1st row)        */
                    1,               /* I - first vector element to write (1 = 1st) */
                    1,               /* I - number of strings to write              */
                    [CString::new("X").unwrap().into_raw()].as_mut_ptr(), /* I - array of pointers to strings            */
                    &mut status, /* IO - error status                           */
                );
                fits_check_status(status)?;

                // POLAA. ffpcle = fits_write_col_flt
                fitsio_sys::ffpcle(
                    uvfits.as_raw(),    /* I - FITS file pointer                       */
                    7,                  /* I - number of column to write (1 = 1st col) */
                    row,                /* I - first row to write (1 = 1st row)        */
                    1,                  /* I - first vector element to write (1 = 1st) */
                    1,                  /* I - number of values to write               */
                    [0.0].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,        /* IO - error status                           */
                );
                fits_check_status(status)?;

                // POL calA
                fitsio_sys::ffpcle(
                    uvfits.as_raw(),    /* I - FITS file pointer                       */
                    8,                  /* I - number of column to write (1 = 1st col) */
                    row,                /* I - first row to write (1 = 1st row)        */
                    1,                  /* I - first vector element to write (1 = 1st) */
                    1,                  /* I - number of values to write               */
                    [0.0].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,        /* IO - error status                           */
                );
                fits_check_status(status)?;

                // POLTYB
                fitsio_sys::ffpcls(
                    uvfits.as_raw(), /* I - FITS file pointer                       */
                    9,               /* I - number of column to write (1 = 1st col) */
                    row,             /* I - first row to write (1 = 1st row)        */
                    1,               /* I - first vector element to write (1 = 1st) */
                    1,               /* I - number of strings to write              */
                    [CString::new("Y").unwrap().into_raw()].as_mut_ptr(), /* I - array of pointers to strings            */
                    &mut status, /* IO - error status                           */
                );
                fits_check_status(status)?;

                // POLAB.
                fitsio_sys::ffpcle(
                    uvfits.as_raw(),     /* I - FITS file pointer                       */
                    10,                  /* I - number of column to write (1 = 1st col) */
                    row,                 /* I - first row to write (1 = 1st row)        */
                    1,                   /* I - first vector element to write (1 = 1st) */
                    1,                   /* I - number of values to write               */
                    [90.0].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,         /* IO - error status                           */
                );
                fits_check_status(status)?;

                // POL calB
                fitsio_sys::ffpcle(
                    uvfits.as_raw(),    /* I - FITS file pointer                       */
                    11,                 /* I - number of column to write (1 = 1st col) */
                    row,                /* I - first row to write (1 = 1st row)        */
                    1,                  /* I - first vector element to write (1 = 1st) */
                    1,                  /* I - number of values to write               */
                    [0.0].as_mut_ptr(), /* I - array of values to write                */
                    &mut status,        /* IO - error status                           */
                );
                fits_check_status(status)?;
            }
        }

        Ok(())
    }

    /// Write a visibility row into the uvfits file.
    ///
    /// `uvfits` must have been opened in write mode and currently have HDU 0
    /// open. The [FitsFile] must be supplied to this function to force the
    /// caller to think about calling this function efficiently; opening the
    /// file for every call would be a problem, and keeping the file open in
    /// [UvfitsWriter] would mean the struct is not thread safe.
    ///
    /// `tile_index1` and `tile_index2` are expected to be zero indexed; they
    /// are made one indexed by this function.
    // TODO: Assumes that all fine channels are written in `vis.` This needs to
    // be updated to add visibilities to an existing uvfits row.
    #[inline]
    pub(crate) fn write_vis(
        &mut self,
        uvfits: &mut FitsFile,
        uvw: &UVW,
        tile_index1: usize,
        tile_index2: usize,
        epoch: &Epoch,
        vis: &[f32],
    ) -> Result<(), UvfitsWriteError> {
        if self.current_num_rows + 1 > self.total_num_rows {
            return Err(UvfitsWriteError::BadRowNum {
                row_num: self.current_num_rows,
                num_rows: self.total_num_rows,
            });
        }

        let mut row = Vec::with_capacity(5 + vis.len());
        row.push((uvw.u / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S) as f32);
        row.push((uvw.v / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S) as f32);
        row.push((uvw.w / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S) as f32);
        row.push(encode_uvfits_baseline(tile_index1 + 1, tile_index2 + 1) as f32);
        let jd_trunc = self.start_epoch.as_jde_utc_days().floor() + 0.5;
        let jd_frac = epoch.as_jde_utc_days() - jd_trunc;
        row.push(jd_frac as f32);
        for &v in vis {
            row.push(v);
        }

        let mut status = 0;
        unsafe {
            fitsio_sys::ffpgpe(
                uvfits.as_raw(),                  /* I - FITS file pointer                      */
                self.current_num_rows as i64 + 1, /* I - group to write(1 = 1st group)          */
                1,                                /* I - first vector element to write(1 = 1st) */
                row.len() as i64,                 /* I - number of values to write              */
                row.as_mut_ptr(),                 /* I - array of values that are written       */
                &mut status,                      /* IO - error status                          */
            );
        }
        fits_check_status(status)?;
        self.current_num_rows += 1;
        Ok(())
    }

    // /// Assumes that `vis_array` has already had `weights` applied; these need
    // /// to be undone locally by this function.
    // // TODO: Assumes that all fine channels are written for all baselines in
    // // `vis_array.`
    // pub(crate) fn write_from_vis(
    //     &mut self,
    //     uvfits: &mut FitsFile,
    //     vis_array: ArrayView2<Jones<f32>>,
    //     weights: ArrayView2<f32>,
    //     uvws: &[UVW],
    //     epoch: &Epoch,
    //     num_fine_chans: usize,
    //     fine_chan_flags: &HashSet<usize>,
    // ) -> Result<(), UvfitsWriteError> {
    //     let num_unflagged_baselines = vis_array.len_of(Axis(0));
    //     let num_unflagged_tiles = num_tiles_from_num_baselines(num_unflagged_baselines);
    //     // Write out all the baselines of the timestep we received.
    //     let mut vis: Vec<f32> = Vec::with_capacity(12 * num_fine_chans);
    //     for (unflagged_bl, uvw) in (0..num_unflagged_baselines).into_iter().zip(uvws.iter()) {
    //         // uvfits expects the tile numbers to be from 1 to the total number
    //         // of tiles sequentially, so don't use the actual unflagged tile
    //         // numbers.
    //         let (tile1, tile2) =
    //             cross_correlation_baseline_to_tiles(num_unflagged_tiles, unflagged_bl);
    //         let mut unflagged_chan_index = 0;
    //         for fine_chan_index in 0..num_fine_chans {
    //             if fine_chan_flags.contains(&fine_chan_index) {
    //                 vis.extend_from_slice(&[0.0; 12])
    //             } else {
    //                 let weight = unsafe { weights.uget((unflagged_bl, unflagged_chan_index)) };
    //                 let jones = (unsafe { vis_array.uget((unflagged_bl, unflagged_chan_index)) })
    //                     .clone()
    //                     // Undo the weight.
    //                     * (1.0 / weight);
    //                 unflagged_chan_index += 1;
    //                 vis.extend_from_slice(&[
    //                     // XX
    //                     jones[0].re,
    //                     jones[0].im,
    //                     *weight,
    //                     // YY
    //                     jones[3].re,
    //                     jones[3].im,
    //                     *weight,
    //                     // XY
    //                     jones[1].re,
    //                     jones[1].im,
    //                     *weight,
    //                     // YX
    //                     jones[2].re,
    //                     jones[2].im,
    //                     *weight,
    //                 ]);
    //             };
    //         }
    //         self.write_vis(uvfits, uvw, tile1, tile2, epoch, &vis)?;
    //         vis.clear();
    //     }
    //     Ok(())
    // }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mwalib::{
        _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
        get_required_fits_key,
    };
    use tempfile::NamedTempFile;

    use float_cmp::{approx_eq, F64Margin};

    use crate::{get_flaggable_timesteps, constants::{COTTER_MWA_LATITUDE_RADIANS, HIFITIME_GPS_FACTOR}};

    #[test]
    // Make a tiny uvfits file. The result has been verified by CASA's
    // "importuvfits" function.
    fn test_new_uvfits_is_sensible() {
        let tmp_uvfits_file = NamedTempFile::new().unwrap();
        let num_timesteps = 1;
        let num_baselines = 3;
        let num_chans = 2;
        let obsid = 1065880128;
        let start_epoch = Epoch::from_tai_seconds(obsid as f64 + 19.0 + HIFITIME_GPS_FACTOR);

        let mut u = UvfitsWriter::new(
            // tmp_uvfits_file.path(),
            Path::new("tests/data/test.uvfits"),
            num_timesteps,
            num_baselines,
            num_chans,
            &start_epoch,
            40e3,
            170e6,
            3,
            &RADec::new_degrees(0.0, 60.0),
            Some("test"),
        )
        .unwrap();

        let mut f = u.open().unwrap();
        for _timestep_index in 0..num_timesteps {
            for baseline_index in 0..num_baselines {
                let (tile1, tile2) = match baseline_index {
                    0 => (0, 1),
                    1 => (0, 2),
                    2 => (1, 2),
                    _ => unreachable!(),
                };

                u.write_vis(
                    &mut f,
                    &UVW::default(),
                    tile1,
                    tile2,
                    &start_epoch,
                    (baseline_index..baseline_index + num_chans)
                        .into_iter()
                        .map(|int| int as f32)
                        .collect::<Vec<_>>()
                        .as_slice(),
                )
                .unwrap();
            }
        }

        let names = ["Tile1", "Tile2", "Tile3"];
        let positions: Vec<XyzGeodetic> = (0..names.len())
            .into_iter()
            .map(|i| XyzGeodetic {
                x: i as f64,
                y: i as f64 * 2.0,
                z: i as f64 * 3.0,
            })
            .collect();
        u.write_uvfits_antenna_table(&names, &positions).unwrap();
    }

    // TODO: dedup this from lib.rs
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

    fn assert_uvfits_primary_hdu_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut left_primary_hdu = fits_open_hdu!(left_fptr, 0).unwrap();
        let mut right_primary_hdu = fits_open_hdu!(right_fptr, 0).unwrap();

        // WONTFIX:
        // "METAVER",
        // "COTVER"
        // "MWAPYVER",
        // "OBJECT",

        let short_string_keys = vec![
            "SIMPLE", "EXTEND", "GROUPS", "PCOUNT", "GCOUNT", "PTYPE1", "PTYPE2", "PTYPE3",
            "PTYPE4", "PTYPE5", "CTYPE2", "CTYPE3", "CTYPE4", "CTYPE5", "CTYPE6", "TELESCOP", 
            "INSTRUME",
            // TODO: "DATE-OBS",
        ];

        let f64_keys = vec![
            "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "NAXIS4", "NAXIS5", "NAXIS6",
            "BSCALE", "PSCAL1", "PZERO1", "PSCAL2", "PZERO2", "PSCAL3", "PZERO3", "PSCAL4",
            "PZERO4", "PSCAL5", "PZERO5", "CRVAL2", "CRPIX2", "CDELT2", "CRVAL3", "CDELT3",
            "CRPIX3", "CRVAL4", "CDELT4", "CRPIX4", "CRVAL5", "CRPIX5", "CDELT5", "CRVAL6",
            "CRPIX6", "CDELT6", "EPOCH", "OBSRA", "OBSDEC",
            // TODO: "FIBRFACT",
        ];

        for key in short_string_keys {
            let left_result: Result<String, _> =
                get_required_fits_key!(left_fptr, &left_primary_hdu, key);
            let right_result: Result<String, _> =
                get_required_fits_key!(right_fptr, &right_primary_hdu, key);
            match (left_result, right_result) {
                (Ok(left_val), Ok(right_val)) => {
                    assert_eq!(
                        left_val, right_val,
                        "mismatch for metafits short string key {}",
                        key,
                    );
                }
                (Err(err), Ok(right_val)) => {
                    panic!(
                        "unable to get left short string key {}. Right val={:?}. err={}",
                        key, right_val, err
                    );
                }
                (.., Err(err)) => {
                    panic!("unable to get right short string key {}. {}", key, err);
                }
            }
        }

        for key in f64_keys {
            let left_result: Result<f64, _> =
                get_required_fits_key!(left_fptr, &left_primary_hdu, key);
            let right_result: Result<f64, _> =
                get_required_fits_key!(right_fptr, &right_primary_hdu, key);
            match (left_result, right_result) {
                (Ok(left_val), Ok(right_val)) => {
                    assert!(
                        approx_eq!(f64, left_val, right_val, F64Margin::default()),
                        "mismatch for metafits short f64 key {}. {} != {} (margin={:?})",
                        key,
                        left_val,
                        right_val,
                        F64Margin::default()
                    );
                }
                (Err(err), Ok(right_val)) => {
                    panic!(
                        "unable to get left short f64 key {}. Right val={}. err={}",
                        key, right_val, err
                    );
                }
                (.., Err(err)) => {
                    panic!("unable to get right short f64 key {}. {}", key, err);
                }
            }
        }
    }

    fn assert_uvfits_ants_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut left_ant_hdu = fits_open_hdu!(left_fptr, 1).unwrap();
        let mut right_ant_hdu = fits_open_hdu!(right_fptr, 1).unwrap();

        panic!("TODO: read ants");
    }

    #[test]
    fn test_uvfits_from_context_matches_cotter_header() {
        let context = get_mwa_ord_context();

        let obsid = context.metafits_context.obs_id;
        let start_epoch = Epoch::from_tai_seconds(obsid as f64 + 19.0 + HIFITIME_GPS_FACTOR);

        let tmp_uvfits_file = NamedTempFile::new().unwrap();
        // let tmp_uvfits_file = Path::new("tests/data/test_header.uvfits");

        let timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let coarse_chan_idxs = context.common_coarse_chan_indices.clone();

        let mut u = UvfitsWriter::from_context(
            tmp_uvfits_file.path(),
            &context,
            &timestep_idxs,
            &coarse_chan_idxs,
            &start_epoch,
        )
        .unwrap();

        let mut f = u.open().unwrap();

        drop(f);

        let cotter_uvfits_path = Path::new("tests/data/1196175296_mwa_ord/1196175296.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_primary_hdu_eq(&mut birli_fptr, &mut cotter_fptr)
    }

    
    #[test]
    fn test_uvfits_vis_from_imgsets_matches_cotter() {

        let context = get_mwa_ord_context();

        let obsid = context.metafits_context.obs_id;
        let start_epoch = Epoch::from_tai_seconds(obsid as f64 + 19.0 + HIFITIME_GPS_FACTOR);

        let tmp_uvfits_file = Path::new("tests/data/test_ants.uvfits");

        let timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let coarse_chan_idxs = context.common_coarse_chan_indices.clone();

        let mut u = UvfitsWriter::from_context(
            tmp_uvfits_file,
            &context,
            &timestep_idxs,
            &coarse_chan_idxs,
            &start_epoch,
        )
        .unwrap();

        let mut f = u.open().unwrap();

    }

    #[test]
    fn test_uvfits_antenna_table_from_context_matches_cotter() {

        let context = get_mwa_ord_context();

        let obsid = context.metafits_context.obs_id;
        let start_epoch = Epoch::from_tai_seconds(obsid as f64 + 19.0 + HIFITIME_GPS_FACTOR);

        let tmp_uvfits_file = Path::new("tests/data/test_ants.uvfits");

        let timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let coarse_chan_idxs = context.common_coarse_chan_indices.clone();

        let mut u = UvfitsWriter::from_context(
            tmp_uvfits_file,
            &context,
            &timestep_idxs,
            &coarse_chan_idxs,
            &start_epoch,
        )
        .unwrap();

        let mut f = u.open().unwrap();

        panic!("TODO: write visibilities!");
        
        u.write_ants_from_context(&context, Some(COTTER_MWA_LATITUDE_RADIANS)).unwrap();

        drop(f);

        let cotter_uvfits_path = Path::new("tests/data/1196175296_mwa_ord/1196175296.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_ants_eq(&mut birli_fptr, &mut cotter_fptr)

        // TODO
        // let tmp_uvfits_file = NamedTempFile::new().unwrap();
        // let num_timesteps = 1;
        // let num_baselines = 3;
        // let num_chans = 2;
        // let obsid = 1065880128;
        // let start_epoch = Epoch::from_tai_seconds(obsid as f64 + 19.0 + HIFITIME_GPS_FACTOR);

        // let mut u = UvfitsWriter::new(
        //     // tmp_uvfits_file.path(),
        //     Path::new("tests/data/test_antenna.uvfits"),
        //     num_timesteps,
        //     num_baselines,
        //     num_chans,
        //     &start_epoch,
        //     40e3,
        //     170e6,
        //     3,
        //     &RADec::new_degrees(0.0, 60.0),
        //     Some("test"),
        // )
        // .unwrap();

        // let mut f = u.open().unwrap();
        // for _timestep_index in 0..num_timesteps {
        //     for baseline_index in 0..num_baselines {
        //         let (tile1, tile2) = match baseline_index {
        //             0 => (0, 1),
        //             1 => (0, 2),
        //             2 => (1, 2),
        //             _ => unreachable!(),
        //         };

        //         u.write_vis(
        //             &mut f,
        //             &UVW::default(),
        //             tile1,
        //             tile2,
        //             &start_epoch,
        //             (baseline_index..baseline_index + num_chans)
        //                 .into_iter()
        //                 .map(|int| int as f32)
        //                 .collect::<Vec<_>>()
        //                 .as_slice(),
        //         )
        //         .unwrap();
        //     }
        // }

        // let names = ["Tile1", "Tile2", "Tile3"];
        // let positions: Vec<XyzGeodetic> = (0..names.len())
        //     .into_iter()
        //     .map(|i| XyzGeodetic {
        //         x: i as f64,
        //         y: i as f64 * 2.0,
        //         z: i as f64 * 3.0,
        //     })
        //     .collect();
        // u.write_uvfits_antenna_table(&names, &positions).unwrap();
    }
}
