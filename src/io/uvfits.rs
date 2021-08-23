// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Module for uvfits file format reading and writing
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

// use std::collections::HashSet;
use std::ffi::CString;
use std::path::Path;

use fitsio::{errors::check_status as fits_check_status, FitsFile};
use indicatif::ProgressStyle;
use itertools::izip;
use log::warn;
use mwa_rust_core::{
    constants::VEL_C, erfa_sys, erfa_sys::ERFA_DJM0, fitsio, fitsio_sys, hifitime::Epoch, mwalib,
    precession::*, time, LatLngHeight, RADec, XyzGeodetic, ENH, UVW,
};
// use ndarray::prelude::*;

use super::error::UvfitsWriteError;
use crate::cxx_aoflagger::ffi::{CxxFlagMask, CxxImageSet};
use cxx::UniquePtr;
use mwalib::{Baseline, CorrelatorContext, MetafitsContext};

/// From a `hifitime` [Epoch], get a formatted date string with the hours,
/// minutes and seconds set to 0.
fn get_truncated_date_string(epoch: Epoch) -> String {
    let (year, month, day, _, _, _, _) = epoch.as_gregorian_utc();
    format!(
        "{year}-{month:02}-{day:02}T00:00:00.0",
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
#[allow(dead_code)]
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
pub struct UvfitsWriter<'a> {
    /// The path to the uvifts file.
    path: &'a Path,

    /// The number of uvfits rows. This is equal to `num_timesteps` *
    /// `num_baselines`.
    total_num_rows: usize,

    /// The number of uvfits rows that have currently been written.
    current_num_rows: usize,

    /// The center frequency of the center fine channel \[Hz\].
    centre_freq: f64,

    /// A `hifitime` [Epoch] struct associated with the first timestep of the
    /// data.
    start_epoch: Epoch,

    /// The RA/Dec where this observation is phased to
    phase_centre: RADec,

    /// Array Position, geocentric
    array_pos: LatLngHeight,
}

impl<'a> UvfitsWriter<'a> {
    /// Create a new uvfits file at the specified filename.
    ///
    /// This will destroy any existing uvfits file at that path
    ///
    /// `num_timesteps`, `num_baselines` and `num_chans` are the number of
    /// timesteps, baselines and channels in this uvfits file respectively.
    ///
    /// `start_epoch` can be calculated from GPS Time using the hifitime library, e.g.
    ///
    /// ```rust
    /// use mwa_rust_core::time;
    /// let first_gps_time = 1196175296.0;
    /// let start_epoch = time::gps_to_epoch(first_gps_time);
    /// ```
    ///
    /// `centre_freq_hz` is the centre frequency of the coarse band that this
    /// uvfits file pertains to [Hz].
    ///
    /// `centre_freq_chan` is the index (from zero) of the centre frequency of
    /// the channels of this uvfits file.
    ///
    /// `phase_centre` is a [`RADec`] of the observation's phase center, used to
    /// populate the `OBSRA` and `OBSDEC` keys.
    ///
    /// `obs_name` an optional name for the object under observation. Used to
    /// populate the `OBJECT` keys.
    ///
    /// # Errors
    ///
    /// Will return an [`UvfitsWriteError`] if:
    /// - there is an existing file at `filename` which cannot be removed.
    /// - a fits operation fails.
    ///
    /// TODO: replace all these args with birli_context
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        filename: &'a Path,
        num_timesteps: usize,
        num_baselines: usize,
        num_chans: usize,
        start_epoch: Epoch,
        fine_chan_width_hz: f64,
        centre_freq_hz: f64,
        centre_freq_chan: usize,
        phase_centre: RADec,
        obs_name: Option<&str>,
        array_pos: Option<LatLngHeight>,
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
        hdu.write_key(&mut u, "DATE-OBS", get_truncated_date_string(start_epoch))?;

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

        // TODO: write obsid?

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

        let array_pos = match array_pos {
            Some(pos) => pos,
            None => {
                warn!("we are using MWA lat / lng / height from mwalib for array position");
                // The results here are slightly different to those given by cotter.
                // This is at least partly due to different constants (the altitude is
                // definitely slightly different), but possibly also because ERFA is
                // more accurate than cotter's "homebrewed" Geodetic2XYZ.
                LatLngHeight::new_mwa()
            }
        };

        Ok(Self {
            path: filename,
            total_num_rows,
            current_num_rows: 0,
            centre_freq: centre_freq_hz,
            start_epoch,
            phase_centre,
            array_pos,
        })
    }

    /// Create a new uvfits file at the specified filename using an [`mwalib::CorrelatorContext`]
    ///
    /// start epoch is determined by the first time provided in `img_timestep_idxs`
    /// which may not necessarily match the obsid.
    ///
    /// # Errors
    ///
    /// See: [`UvfitsWriter::new`]
    pub fn from_mwalib(
        filename: &'a Path,
        context: &CorrelatorContext,
        img_timestep_idxs: &[usize],
        img_coarse_chan_idxs: &[usize],
        img_baseline_idxs: &[usize],
        array_pos: Option<LatLngHeight>,
    ) -> Result<Self, UvfitsWriteError> {
        let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let num_img_coarse_chans = img_coarse_chan_idxs.len();
        let num_img_chans = num_fine_chans_per_coarse * num_img_coarse_chans;
        let first_gps_time = context.timesteps[img_timestep_idxs[0]].gps_time_ms as f64 / 1000.0;
        let start_epoch = time::gps_to_epoch(first_gps_time);
        let phase_centre = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
        let centre_freq_chan = num_img_chans / 2;
        Self::new(
            filename,
            img_timestep_idxs.len(),
            img_baseline_idxs.len(),
            num_img_chans,
            start_epoch,
            context.metafits_context.corr_fine_chan_width_hz as f64,
            context.metafits_context.centre_freq_hz as f64,
            centre_freq_chan,
            phase_centre,
            Some(context.metafits_context.obs_name.as_str()),
            array_pos,
        )
    }

    /// Opens the associated uvfits file in edit mode, returning the [FitsFile]
    /// struct.
    ///
    /// Closing the file can be achieved by `drop`ing all references to the
    /// [`FitsFile`]
    ///
    /// # Errors
    ///
    /// Will return an [`fitsio::errors::Error`] if the file can't be edited.
    pub fn open(&self) -> Result<FitsFile, fitsio::errors::Error> {
        let mut f = FitsFile::edit(&self.path)?;
        // Ensure HDU 0 is opened.
        f.hdu(0)?;
        Ok(f)
    }

    /// Write the antenna table to a uvfits file. Assumes that the array
    /// location is MWA.
    ///
    /// `positions` are the [XyzGeodetic] coordinates
    /// of the MWA tiles. These positions need to have the MWA's "centre" XYZ
    /// coordinates subtracted to make them local XYZ.
    ///
    /// `Self` must have only have a single HDU when this function is called
    /// (true when using methods only provided by `Self`).
    ///
    /// Derived from cotter.
    ///
    /// # Errors
    ///
    /// Will return an [`UvfitsWriteError`] if a fits operation fails.
    pub fn write_uvfits_antenna_table<T: AsRef<str>>(
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

        let array_xyz = self.array_pos.to_geocentric_wgs84()?;

        hdu.write_key(&mut uvfits, "ARRAYX", array_xyz.x)?;
        hdu.write_key(&mut uvfits, "ARRAYY", array_xyz.y)?;
        hdu.write_key(&mut uvfits, "ARRAYZ", array_xyz.z)?;

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

    /// Write the antenna table to a uvfits file using the provided
    /// [`mwalib::CorrelatorContext`]
    ///
    /// `self` must have only have a single HDU when this function is called
    /// (true when using methods only provided by `Self`).
    ///
    /// `latitude_rad` Optionally override the latitude of the array [Radians]
    ///
    /// # Errors
    ///
    /// See: [`UvfitsWriter::write_uvfits_antenna_table`]
    ///
    pub fn write_ants_from_mwalib(self, context: &MetafitsContext) -> Result<(), UvfitsWriteError> {
        let (antenna_names, positions): (Vec<String>, Vec<XyzGeodetic>) = context
            .antennas
            .iter()
            .map(|antenna| {
                let position_enh = ENH {
                    e: antenna.east_m,
                    n: antenna.north_m,
                    h: antenna.height_m,
                };
                let position = position_enh.to_xyz(self.array_pos.latitude_rad);
                (antenna.tile_name.to_owned(), position)
            })
            .unzip();
        self.write_uvfits_antenna_table(&antenna_names, &positions)
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
    ///
    /// # Errors
    ///
    /// Will return an [`UvfitsWriteError`] if a fits operation fails.
    ///
    /// TODO: Assumes that all fine channels are written in `vis.` This needs to
    /// be updated to add visibilities to an existing uvfits row.
    ///
    /// TODO: replace all these args with birli_context
    #[allow(clippy::too_many_arguments)]
    #[inline]
    pub fn write_vis(
        &mut self,
        uvfits: &mut FitsFile,
        uvw: UVW,
        tile_index1: usize,
        tile_index2: usize,
        epoch: Epoch,
        vis: &[f32],
    ) -> Result<(), UvfitsWriteError> {
        if self.current_num_rows + 1 > self.total_num_rows {
            return Err(UvfitsWriteError::BadRowNum {
                row_num: self.current_num_rows,
                num_rows: self.total_num_rows,
            });
        }

        let mut row = Vec::with_capacity(5 + vis.len());
        row.push((uvw.u / VEL_C) as f32);
        row.push((uvw.v / VEL_C) as f32);
        row.push((uvw.w / VEL_C) as f32);
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

    /// Write visibilty rows into the uvfits file from the provided
    /// [`mwalib::CorrelatorContext`].
    ///
    /// `uvfits` must have been opened in write mode and currently have HDU 0
    /// open. The [FitsFile] must be supplied to this function to force the
    /// caller to think about calling this function efficiently; opening the
    /// file for every call would be a problem, and keeping the file open in
    /// [UvfitsWriter] would mean the struct is not thread safe.
    ///
    /// `baseline_idxs` the baseline indices (according to mwalib)
    /// which should be written to the file
    ///
    /// `baseline_imgsets` a [`CxxImageSet`] for each baseline in `baseline_idxs`
    ///
    /// `baseline_flagmasks` a [`CxxFlagMask`] for each baseline in `baseline_idxs`
    ///
    /// `img_timestep_idxs` the timestep indices (according to mwalib)
    /// which are used in all of the images in `baseline_imgsets` and all of the
    /// flagmasks in `baseline_flagmasks`
    ///
    /// `img_coarse_chan_idxs` the coarse channel indices (according to mwalib)
    /// which are used in all of the images in `baseline_imgsets` and all of the
    /// flagmasks in `baseline_flagmasks`
    ///
    /// TODO: handle averaging
    ///
    /// # Errors
    ///
    /// Will return an [`UvfitsWriteError`] if a fits operation fails.
    ///
    ///
    /// TODO: replace all these args with birli_context
    #[allow(clippy::too_many_arguments)]
    pub fn write_baseline_imgset_flagmasks(
        &mut self,
        uvfits: &mut FitsFile,
        context: &CorrelatorContext,
        baseline_idxs: &[usize],
        baseline_imgsets: &[UniquePtr<CxxImageSet>],
        baseline_flagmasks: &[UniquePtr<CxxFlagMask>],
        img_timestep_idxs: &[usize],
        img_coarse_chan_idxs: &[usize],
    ) -> Result<(), UvfitsWriteError> {
        let num_baselines = baseline_idxs.len();
        assert_eq!(num_baselines, baseline_imgsets.len());
        assert_eq!(num_baselines, baseline_flagmasks.len());

        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);

        let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let num_img_coarse_chans = img_coarse_chan_idxs.len();
        let num_img_chans = num_fine_chans_per_coarse * num_img_coarse_chans;

        let num_img_timesteps = img_timestep_idxs.len();
        assert_eq!(self.total_num_rows, num_img_timesteps * num_baselines);

        // Weights are normalized so that default res of 10 kHz, 1s has weight of "1" per sample

        // TODO: deal with weight factor when doing averaging.
        // TODO: deal with passband gains
        let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;
        let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz as f64;
        let weight_factor = fine_chan_width_hz * integration_time_s / 10000.0;

        // MWA/CASA/AOFlagger visibility order is XX,XY,YX,YY
        // UVFits visibility order is XX,YY,XY,YX
        let pol_order = vec![0, 3, 1, 2];
        // Real and imaginary component offsets within image width
        // let re_im = vec![0, 1];

        let num_blts = img_timestep_idxs.len() * baseline_idxs.len();

        // Create a progress bar to show the writing status
        let write_progress = indicatif::ProgressBar::new(num_blts as u64);
        write_progress.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .progress_chars("=> "),
        );
        write_progress.set_message("write uv vis");

        for (img_timestep_idx, &timestep_idx) in img_timestep_idxs.iter().enumerate() {
            let gps_time_s = context.timesteps[timestep_idx].gps_time_ms as f64 / 1000.0;
            let epoch = time::gps_to_epoch(gps_time_s + integration_time_s / 2.0);

            let prec_info = precess_time(
                self.phase_centre,
                epoch,
                self.array_pos.longitude_rad,
                self.array_pos.latitude_rad,
            );

            let tiles_xyz_precessed = prec_info.precess_xyz_parallel(&tiles_xyz_geod);

            for (&baseline_idx, imgset, flagmask) in
                izip!(baseline_idxs, baseline_imgsets, baseline_flagmasks)
            {
                let baseline: &Baseline = &context.metafits_context.baselines[baseline_idx];
                let imgset: &UniquePtr<CxxImageSet> = imgset;
                let img_stride: usize = imgset.HorizontalStride();

                let flagmask: &UniquePtr<CxxFlagMask> = flagmask;
                let flag_buffer: &[bool] = flagmask.Buffer();
                let flag_stride: usize = flagmask.HorizontalStride();

                let ant1_idx = baseline.ant1_index;
                let ant2_idx = baseline.ant2_index;

                let baseline_xyz_precessed =
                    tiles_xyz_precessed[ant1_idx] - tiles_xyz_precessed[ant2_idx];
                let uvw = UVW::from_xyz(baseline_xyz_precessed, prec_info.hadec_j2000);

                // TODO: this is extremely inefficient.

                let vis: Vec<f32> = (0..num_img_chans)
                    .flat_map(|chan_idx| {
                        pol_order
                            .iter()
                            .flat_map(|img_pol_idx| {
                                let img_buffer_re: &[f32] = imgset.ImageBuffer(img_pol_idx * 2);
                                let img_buffer_im: &[f32] = imgset.ImageBuffer(img_pol_idx * 2 + 1);
                                let vis_re =
                                    img_buffer_re[chan_idx * img_stride + img_timestep_idx];
                                let vis_im =
                                    img_buffer_im[chan_idx * img_stride + img_timestep_idx];
                                let flag = flag_buffer[chan_idx * flag_stride + img_timestep_idx];
                                // TODO: weight from flags
                                let mut weight = weight_factor as f32;
                                if flag {
                                    weight *= -1.0;
                                }
                                vec![vis_re, vis_im, weight]
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect();

                // TODO: calculate weights
                let result = self.write_vis(uvfits, uvw, ant1_idx, ant2_idx, epoch, &vis.clone());
                if let Err(err) = result {
                    return Err(err);
                };

                write_progress.inc(1);
            }
        }

        write_progress.finish();

        Ok(())
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use mwa_rust_core::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib, time,
    };
    use mwalib::{
        _get_fits_col, _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
        get_fits_col, get_required_fits_key,
    };
    use tempfile::NamedTempFile;

    use float_cmp::{approx_eq, F32Margin, F64Margin};

    use crate::{
        context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets_existing, get_antenna_flags,
        get_flaggable_timesteps, init_baseline_flagmasks,
    };

    use fitsio::{
        errors::check_status as fits_check_status,
        hdu::{FitsHdu, HduInfo},
    };

    #[test]
    // Make a tiny uvfits file. The result has been verified by CASA's
    // "importuvfits" function.
    fn test_new_uvfits_is_sensible() {
        let tmp_uvfits_file = NamedTempFile::new().unwrap();
        let num_timesteps = 1;
        let num_baselines = 3;
        let num_chans = 2;
        let obsid = 1065880128;
        let start_epoch = time::gps_to_epoch(obsid as f64);

        let mut u = UvfitsWriter::new(
            tmp_uvfits_file.path(),
            num_timesteps,
            num_baselines,
            num_chans,
            start_epoch,
            40e3,
            170e6,
            3,
            RADec::new_degrees(0.0, 60.0),
            Some("test"),
            None,
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
                    UVW::default(),
                    tile1,
                    tile2,
                    start_epoch,
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

    macro_rules! assert_short_string_keys_eq {
        ($keys:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr) => {
            for key in $keys {
                match (
                    get_required_fits_key!($left_fptr, &$left_hdu, key),
                    get_required_fits_key!($right_fptr, &$right_hdu, key),
                ) {
                    (Ok::<String, _>(left_val), Ok::<String, _>(right_val)) => {
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
        };
    }

    macro_rules! assert_f64_string_keys_eq {
        ($keys:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr) => {
            for key in $keys {
                match (
                    get_required_fits_key!($left_fptr, &$left_hdu, key),
                    get_required_fits_key!($right_fptr, &$right_hdu, key),
                ) {
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
        };
    }

    macro_rules! assert_table_column_descriptions_match {
        ( $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr ) => {
            match (&$left_hdu.info, &$right_hdu.info) {
                (
                    HduInfo::TableInfo {
                        column_descriptions: left_columns,
                        ..
                    },
                    HduInfo::TableInfo {
                        column_descriptions: right_columns,
                        ..
                    },
                ) => {
                    for (col_idx, (left_col, right_col)) in
                        izip!(left_columns, right_columns).enumerate()
                    {
                        assert_eq!(
                            left_col, right_col,
                            "column description at index {} does not match",
                            col_idx
                        );
                    }
                    left_columns
                        .iter()
                        .map(|col| col.name.clone())
                        .collect::<Vec<String>>()
                }
                _ => {
                    panic!("could not read left ant HDU as a table")
                }
            };
        };
    }

    macro_rules! assert_table_string_column_values_match {
        ( $col_names:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr ) => {
            for col_name in $col_names {
                let left_col: Vec<String> =
                    get_fits_col!($left_fptr, &$left_hdu, col_name).unwrap();
                let right_col: Vec<String> =
                    get_fits_col!($right_fptr, &$right_hdu, col_name).unwrap();
                assert_eq!(
                    left_col.len(),
                    right_col.len(),
                    "tables not the same length."
                );
                for (row_idx, (left_cell, right_cell)) in izip!(left_col, right_col).enumerate() {
                    assert_eq!(
                        left_cell, right_cell,
                        "cells don't match in column {}, row {}",
                        col_name, row_idx
                    );
                }
            }
        };
    }

    macro_rules! assert_table_f64_column_values_match {
        ( $col_names:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr ) => {
            for col_name in $col_names {
                let left_col: Vec<f64> = get_fits_col!($left_fptr, &$left_hdu, col_name).unwrap();
                let right_col: Vec<f64> =
                    get_fits_col!($right_fptr, &$right_hdu, col_name).unwrap();
                assert_eq!(
                    left_col.len(),
                    right_col.len(),
                    "tables not the same length."
                );
                for (row_idx, (left_cell, right_cell)) in izip!(left_col, right_col).enumerate() {
                    assert!(
                        approx_eq!(f64, left_cell, right_cell, F64Margin::default()),
                        "cells don't match in column {}, row {}. {} != {}",
                        col_name,
                        row_idx,
                        left_cell,
                        right_cell
                    );
                }
            }
        };
    }

    macro_rules! assert_table_vector_f64_column_values_match {
        ( $col_descriptions:expr, $col_info:expr, $row_len:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr ) => {
            for (col_name, len) in $col_info
                .into_iter()
                .map(|(_str, len)| (_str.to_string(), len))
            {
                let col_num = $col_descriptions
                    .iter()
                    .position(|description| description == &col_name)
                    .expect(format!("could not find a column with the name {}", col_name).as_str());

                let mut status = 0;
                let mut left_cell_vector: Vec<f64> = vec![0.0; len as usize];
                let mut right_cell_vector: Vec<f64> = vec![0.0; len as usize];

                for row_idx in 0..$row_len {
                    unsafe {
                        fitsio_sys::ffgcvd(
                            $left_fptr.as_raw(),
                            (col_num + 1) as _,
                            (row_idx + 1) as _,
                            1,
                            len,
                            0 as _,
                            left_cell_vector.as_mut_ptr() as _,
                            &mut 0,
                            &mut status,
                        );
                        fits_check_status(status).unwrap();
                        fitsio_sys::ffgcvd(
                            $right_fptr.as_raw(),
                            (col_num + 1) as _,
                            (row_idx + 1) as _,
                            1,
                            len,
                            0 as _,
                            right_cell_vector.as_mut_ptr() as _,
                            &mut 0,
                            &mut status,
                        );
                        fits_check_status(status).unwrap();
                    }
                    for (cell_idx, (&left_cell, &right_cell)) in
                        izip!(&left_cell_vector, &right_cell_vector).enumerate()
                    {
                        assert!(
                            approx_eq!(
                                f64,
                                left_cell,
                                right_cell,
                                F64Margin::default().epsilon(1e-12)
                            ),
                            "cells don't match in column {}, row {}, cell index {}. {} != {}",
                            col_name,
                            row_idx,
                            cell_idx,
                            left_cell,
                            right_cell
                        );
                    }
                }
            }
        };
    }

    fn assert_uvfits_primary_header_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut _left_primary_hdu = fits_open_hdu!(left_fptr, 0).unwrap();
        let mut _right_primary_hdu = fits_open_hdu!(right_fptr, 0).unwrap();

        assert_short_string_keys_eq!(
            vec![
                "SIMPLE", "EXTEND", "GROUPS", "PCOUNT", "GCOUNT", "PTYPE1", "PTYPE2", "PTYPE3",
                "PTYPE4", "PTYPE5", "CTYPE2", "CTYPE3", "CTYPE4", "CTYPE5", "CTYPE6", "TELESCOP",
                "INSTRUME",
                "DATE-OBS",
                // WONTFIX:
                // "METAVER",
                // "COTVER"
                // "MWAPYVER",
                // "OBJECT",
            ],
            left_fptr,
            _left_primary_hdu,
            right_fptr,
            _right_primary_hdu
        );

        assert_f64_string_keys_eq!(
            vec![
                "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "NAXIS4", "NAXIS5", "NAXIS6",
                "BSCALE", "PSCAL1", "PZERO1", "PSCAL2", "PZERO2", "PSCAL3", "PZERO3", "PSCAL4",
                "PZERO4", "PSCAL5", "PZERO5", "CRVAL2", "CRPIX2", "CDELT2", "CRVAL3", "CDELT3",
                "CRPIX3", "CRVAL4", "CDELT4", "CRPIX4", "CRVAL5", "CRPIX5", "CDELT5", "CRVAL6",
                "CRPIX6", "CDELT6", "EPOCH", "OBSRA",
                "OBSDEC",
                // TODO: "FIBRFACT", (v-factor, cable type?)
            ],
            left_fptr,
            _left_primary_hdu,
            right_fptr,
            _right_primary_hdu
        );
    }

    fn assert_uvfits_ant_header_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut _left_ant_hdu = fits_open_hdu!(left_fptr, 1).unwrap();
        let mut _right_ant_hdu = fits_open_hdu!(right_fptr, 1).unwrap();

        assert_short_string_keys_eq!(
            vec![
                "XTENSION", "TTYPE1", "TFORM1", "TTYPE2", "TFORM2", "TUNIT2", "TTYPE3", "TFORM3",
                "TTYPE4", "TFORM4", "TTYPE5", "TFORM5", "TUNIT5", "TTYPE6", "TFORM6", "TTYPE7",
                "TFORM7", "TUNIT7", "TTYPE8", "TFORM8", "TTYPE9", "TFORM9", "TTYPE10", "TFORM10",
                "TUNIT10", "TTYPE11", "TFORM11", "EXTNAME", "TIMSYS", "ARRNAM",
            ],
            left_fptr,
            _left_ant_hdu,
            right_fptr,
            _right_ant_hdu
        );

        assert_f64_string_keys_eq!(
            vec![
                "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS",
                // TODO: "ARRAYX", "ARRAYY", "ARRAYZ",
                "FREQ",   // TODO: "GSTIA0",
                "DEGPDY", // TODO: "RDATE",
                "POLARX", "POLARY", "UT1UTC", "DATUTC", "NUMORB", "NOPCAL", "FREQID", "IATUTC",
            ],
            left_fptr,
            _left_ant_hdu,
            right_fptr,
            _right_ant_hdu
        );
    }

    pub fn get_group_column_description(
        fptr: &mut FitsFile,
        hdu: &FitsHdu,
    ) -> Result<Vec<String>, crate::mwalib::FitsError> {
        let pcount: usize = get_required_fits_key!(fptr, hdu, "PCOUNT")?;
        let mut result = Vec::with_capacity(pcount);
        for p in 0..pcount {
            let ptype: String =
                get_required_fits_key!(fptr, hdu, format!("PTYPE{}", p + 1).as_str())?;
            result.push(ptype);
        }
        Ok(result)
    }

    macro_rules! assert_group_column_descriptions_match {
        ( $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr ) => {
            match (
                get_group_column_description($left_fptr, $left_hdu),
                get_group_column_description($left_fptr, $left_hdu),
            ) {
                (Ok(left_columns), Ok(right_columns)) => {
                    for (col_idx, (left_col, right_col)) in
                        izip!(&left_columns, &right_columns).enumerate()
                    {
                        assert_eq!(
                            left_col, right_col,
                            "column description at index {} does not match",
                            col_idx
                        );
                    }
                    left_columns
                }
                _ => {
                    panic!("could not read HDUs as group table")
                }
            };
        };
    }

    fn assert_uvfits_vis_table_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut _left_vis_hdu = fits_open_hdu!(left_fptr, 0).unwrap();
        let mut _right_vis_hdu = fits_open_hdu!(right_fptr, 0).unwrap();

        let column_info = assert_group_column_descriptions_match!(
            left_fptr,
            &_left_vis_hdu,
            right_fptr,
            &_right_vis_hdu
        );

        let left_length: usize =
            get_required_fits_key!(left_fptr, &_left_vis_hdu, "GCOUNT").unwrap();
        let right_length: usize =
            get_required_fits_key!(right_fptr, &_right_vis_hdu, "GCOUNT").unwrap();

        assert_eq!(left_length, right_length, "group lengths don't match");
        let group_len = left_length;

        let pcount = get_required_fits_key!(left_fptr, &_left_vis_hdu, "PCOUNT").unwrap();
        let floats_per_pol: usize =
            get_required_fits_key!(left_fptr, &_left_vis_hdu, "NAXIS2").unwrap();
        let num_pols: usize = get_required_fits_key!(left_fptr, &_left_vis_hdu, "NAXIS3").unwrap();
        let num_fine_freq_chans: usize =
            get_required_fits_key!(left_fptr, &_left_vis_hdu, "NAXIS4").unwrap();

        let mut left_group_params: Vec<f32> = vec![0.0; pcount];
        let mut right_group_params: Vec<f32> = vec![0.0; pcount];
        let mut left_vis: Vec<f32> = vec![0.0; num_fine_freq_chans * num_pols * floats_per_pol];
        let mut right_vis: Vec<f32> = vec![0.0; num_fine_freq_chans * num_pols * floats_per_pol];
        let mut status = 0;

        let baseline_col_num = column_info
            .iter()
            .position(|name| name == &String::from("BASELINE"))
            .unwrap();

        for row_idx in 0..group_len {
            unsafe {
                // ffggpe = fits_read_grppar_flt
                fitsio_sys::ffggpe(
                    left_fptr.as_raw(),             /* I - FITS file pointer                       */
                    1 + row_idx as i64, /* I - group to read (1 = 1st group)           */
                    1,                  /* I - first vector element to read (1 = 1st)  */
                    pcount as i64,      /* I - number of values to read                */
                    left_group_params.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut status, /* IO - error status                           */
                );
                fits_check_status(status).unwrap();
                // ffggpe = fits_read_grppar_flt
                fitsio_sys::ffggpe(
                    right_fptr.as_raw(),             /* I - FITS file pointer                       */
                    1 + row_idx as i64, /* I - group to read (1 = 1st group)           */
                    1,                  /* I - first vector element to read (1 = 1st)  */
                    pcount as i64,      /* I - number of values to read                */
                    right_group_params.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut status, /* IO - error status                           */
                );
                fits_check_status(status).unwrap();
            }

            for (param_name, left_group_param, right_group_param) in
                izip!(&column_info, &left_group_params, &right_group_params)
            {
                assert!(
                    approx_eq!(
                        f32,
                        *left_group_param,
                        *right_group_param,
                        F32Margin::default()
                    ),
                    "cells don't match in param {}, row {}. {} != {}",
                    param_name,
                    row_idx,
                    left_group_param,
                    right_group_param
                );
            }

            // Don't compare autocorrelations because they're broken.
            let (ant1, ant2) = decode_uvfits_baseline(left_group_params[baseline_col_num] as usize);
            if ant1 == ant2 {
                continue;
            }

            unsafe {
                // ffgpve = fits_read_img_flt
                fitsio_sys::ffgpve(
                    left_fptr.as_raw(),    /* I - FITS file pointer                       */
                    1 + row_idx as i64,    /* I - group to read (1 = 1st group)           */
                    1,                     /* I - first vector element to read (1 = 1st)  */
                    left_vis.len() as i64, /* I - number of values to read                */
                    0.0,                   /* I - value for undefined pixels              */
                    left_vis.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut 0,                /* O - set to 1 if any values are null; else 0 */
                    &mut status,           /* IO - error status                           */
                );
                fits_check_status(status).unwrap();

                // ffgpve = fits_read_img_flt
                fitsio_sys::ffgpve(
                    right_fptr.as_raw(),    /* I - FITS file pointer                       */
                    1 + row_idx as i64,     /* I - group to read (1 = 1st group)           */
                    1,                      /* I - first vector element to read (1 = 1st)  */
                    right_vis.len() as i64, /* I - number of values to read                */
                    0.0,                    /* I - value for undefined pixels              */
                    right_vis.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut 0,                 /* O - set to 1 if any values are null; else 0 */
                    &mut status,            /* IO - error status                           */
                );
                fits_check_status(status).unwrap();
            }

            for (vis_idx, (left_val, right_val)) in izip!(&left_vis, &right_vis).enumerate() {
                assert!(
                    approx_eq!(f32, *left_val, *right_val, F32Margin::default()),
                    "cells don't match in row {}, vis index {}. {:?} != {:?}",
                    row_idx,
                    vis_idx,
                    &left_vis,
                    &right_vis
                );
            }
        }
    }

    fn assert_uvfits_ant_table_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut _left_ant_hdu = fits_open_hdu!(left_fptr, 1).unwrap();
        let mut _right_ant_hdu = fits_open_hdu!(right_fptr, 1).unwrap();

        let column_info = assert_table_column_descriptions_match!(
            left_fptr,
            &_left_ant_hdu,
            right_fptr,
            &_right_ant_hdu
        );

        let first_col_name = &column_info[0];
        let left_length = {
            let left_annames: Vec<String> =
                get_fits_col!(left_fptr, &_left_ant_hdu, first_col_name.as_str()).unwrap();
            left_annames.len()
        };
        let right_length = {
            let right_annames: Vec<String> =
                get_fits_col!(right_fptr, &_right_ant_hdu, first_col_name.as_str()).unwrap();
            right_annames.len()
        };
        assert_eq!(left_length, right_length, "column lengths don't match");
        let row_len = left_length;

        assert_table_string_column_values_match!(
            &["ANNAME", "POLTYA", "POLTYB"],
            left_fptr,
            &_left_ant_hdu,
            right_fptr,
            &_right_ant_hdu
        );

        assert_table_f64_column_values_match!(
            &["NOSTA", "MNTSTA", "STAXOF", "POLAA", "POLAB"],
            left_fptr,
            &_left_ant_hdu,
            right_fptr,
            &_right_ant_hdu
        );

        assert_table_vector_f64_column_values_match!(
            column_info,
            vec![("STABXYZ", 3), ("POLCALA", 3), ("POLCALB", 3)],
            row_len,
            left_fptr,
            &_left_ant_hdu,
            right_fptr,
            &_right_ant_hdu
        );
    }

    #[test]
    fn uvfits_from_mwalib_matches_cotter_header() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();
        // let tmp_uvfits_file = Path::new("tests/data/test_header.uvfits");

        let timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let coarse_chan_idxs = context.common_coarse_chan_indices.clone();

        let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &timestep_idxs,
            &coarse_chan_idxs,
            &baseline_idxs,
            array_pos,
        )
        .unwrap();

        let f = u.open().unwrap();

        drop(f);

        let cotter_uvfits_path = Path::new("tests/data/1196175296_mwa_ord/1196175296.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_primary_header_eq(&mut birli_fptr, &mut cotter_fptr)
    }

    #[test]
    fn uvfits_tables_from_mwalib_matches_cotter() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let img_coarse_chan_idxs = context.common_coarse_chan_indices.clone();

        let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &img_timestep_idxs,
            &img_coarse_chan_idxs,
            &baseline_idxs,
            array_pos,
        )
        .unwrap();

        let mut f = u.open().unwrap();

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(get_antenna_flags(&context)),
        );

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        flag_imgsets_existing(
            &aoflagger,
            strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        u.write_baseline_imgset_flagmasks(
            &mut f,
            &context,
            &baseline_idxs,
            &baseline_imgsets,
            &baseline_flagmasks,
            &img_timestep_idxs,
            &img_coarse_chan_idxs,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        drop(f);

        let cotter_uvfits_path = Path::new("tests/data/1196175296_mwa_ord/1196175296.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_primary_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_vis_table_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_table_eq(&mut birli_fptr, &mut cotter_fptr);
    }
}
