// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Module for uvfits file format reading and writing
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

use std::{
    ffi::CString,
    ops::Range,
    path::{Path, PathBuf},
};

use fitsio::{errors::check_status as fits_check_status, FitsFile};
use indicatif::{ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::{trace, warn};
use marlu::{
    average_chunk_for_pols_f64,
    constants::VEL_C,
    erfa_sys,
    erfa_sys::ERFA_DJM0,
    fitsio, fitsio_sys,
    hifitime::Epoch,
    mwalib::{CorrelatorContext, MetafitsContext},
    precession::*,
    time::{self, gps_millis_to_epoch},
    Jones, LatLngHeight, RADec, XyzGeodetic, ENH, UVW,
};

use crate::{
    flags::{expand_flag_array, flag_to_weight_array, get_weight_factor},
    marlu::num_complex::Complex,
    ndarray::{s, Array, Array3, ArrayView3, ArrayView4, Axis},
};

use super::{
    error::{IOError, UvfitsWriteError},
    WriteableVis,
};

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
pub fn encode_uvfits_baseline(ant1: usize, ant2: usize) -> usize {
    if ant2 > 255 {
        ant1 * 2048 + ant2 + 65_536
    } else {
        ant1 * 256 + ant2
    }
}

/// Decode a uvfits baseline into the antennas that formed it. Antenna indices
/// start at 1.
#[allow(dead_code)]
pub fn decode_uvfits_baseline(bl: usize) -> (usize, usize) {
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
/// Note: only a single contiguous spectral window is supported.
///
/// TODO: make writer and reader a single class?
pub struct UvfitsWriter {
    /// The path to the uvifts file.
    path: PathBuf,

    /// The number of uvfits rows. This is equal to `num_timesteps` *
    /// `num_baselines`.
    total_num_rows: usize,

    /// The number of uvfits rows that have currently been written.
    current_num_rows: usize,

    /// The center frequency of the center fine channel of the spectral
    /// window being written to this file. \[Hz\]
    ///
    /// This is used in both the reference frequency (`FREQ`) in the antenna HDU,
    /// and the center reference value of the frequency axis (`CRVAL4`) in the
    /// visibility hdu.
    centre_freq: f64,

    /// A `hifitime` [Epoch] struct associated with the first timestep of the
    /// data.
    start_epoch: Epoch,

    /// The RA/Dec where this observation is phased to
    phase_centre: RADec,

    /// Array Position [Latitude (radians), Longitude (radians), Height (m)]
    array_pos: LatLngHeight,
}

impl UvfitsWriter {
    /// Create a new uvfits file at the specified path.
    ///
    /// This will destroy any existing uvfits file at that path.
    ///
    /// If you have a [`mwalib::CorrelatorContext`], then it would be more
    /// convenient to use the `from_mwalib` method.
    ///
    /// `num_timesteps`, `num_baselines` and `num_chans` are the number of
    /// timesteps, baselines and channels in this uvfits file respectively. This
    /// is counted after averaging.
    ///
    /// `start_epoch` is a [`hifitime::Epoch`] at the start of the first scan,
    /// and can be calculated from GPS Time using the hifitime library, e.g.
    ///
    /// ```rust
    /// use marlu::time;
    /// let first_gps_time = 1196175296.0;
    /// let start_epoch = time::gps_to_epoch(first_gps_time);
    /// ```
    ///
    /// `centre_freq_hz` is center frequency of the center fine channel of the
    /// spectral window being written to this file. \[Hz\]
    ///
    /// `centre_freq_chan` is the index (from zero) of the center frequency of
    /// the center fine channel of the spectral] window being written to this
    /// file.
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
    /// - there is an existing file at `path` which cannot be removed.
    /// - a fits operation fails.
    ///
    /// TODO: replace all these args with birli_context
    #[allow(clippy::too_many_arguments)]
    pub fn new<T: AsRef<Path>>(
        path: T,
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
        if path.as_ref().exists() {
            trace!("file {:?} exists, deleting", &path.as_ref());
            std::fs::remove_file(&path)?;
        }

        // Create a new fits file.
        let mut status = 0;
        let c_path = CString::new(path.as_ref().to_str().unwrap())?;
        let mut fptr = std::ptr::null_mut();
        trace!(
            "initialising fits file with fitsio_sys ({:?})",
            &path.as_ref()
        );
        unsafe {
            fitsio_sys::ffinit(
                &mut fptr as *mut *mut _, /* O - FITS file pointer                   */
                c_path.as_ptr(),          /* I - name of file to create              */
                &mut status,              /* IO - error status                       */
            );
        }
        fits_check_status(status)?;

        // Initialise the group header. Copied from cotter. -32 means FLOAT_IMG.
        let naxis = 6;
        let mut naxes = [0, 3, 4, num_chans as i64, 1, 1];
        let num_group_params = 5;
        let total_num_rows = num_timesteps * num_baselines;
        assert!(
            total_num_rows > 0,
            "num_timesteps * num_baselines must be > 0"
        );
        trace!("setting group params in fits file ({:?})", &path.as_ref());
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
        trace!("closing fits file ({:?})", &path.as_ref());
        unsafe {
            fitsio_sys::ffclos(fptr, &mut status);
        }
        fits_check_status(status)?;

        // Open the fits file with rust-fitsio.
        trace!("opening  ({:?})", &path.as_ref());
        let mut u = FitsFile::edit(&path)?;
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
            path: path.as_ref().to_path_buf(),
            total_num_rows,
            current_num_rows: 0,
            centre_freq: centre_freq_hz,
            start_epoch,
            phase_centre,
            array_pos,
        })
    }

    /// Create a new uvfits file at the specified path using an [`mwalib::CorrelatorContext`]
    ///
    /// # Details
    ///
    /// start epoch is determined by `mwalib_timestep_range.start` which may not necessarily match
    /// the obsid.
    ///
    /// # Errors
    ///
    /// See: [`UvfitsWriter::new`]
    ///
    /// TODO: fix too_many_arguments
    #[allow(clippy::too_many_arguments)]
    pub fn from_mwalib<T: AsRef<Path>>(
        path: T,
        context: &CorrelatorContext,
        mwalib_timestep_range: &Range<usize>,
        mwalib_coarse_chan_range: &Range<usize>,
        mwalib_baseline_idxs: &[usize],
        array_pos: Option<LatLngHeight>,
        phase_centre: Option<RADec>,
        avg_time: usize,
        avg_freq: usize,
    ) -> Result<Self, UvfitsWriteError> {
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let first_gps_time_s =
            context.timesteps[mwalib_timestep_range.start].gps_time_ms as f64 / 1000.0;
        let start_epoch = time::gps_to_epoch(first_gps_time_s);
        let phase_centre = phase_centre
            .unwrap_or_else(|| RADec::from_mwalib_phase_or_pointing(&context.metafits_context));

        // TODO: move these calculations into Marlu context.
        let pre_avg_start_chan_idx = mwalib_coarse_chan_range.start * fine_chans_per_coarse;
        let pre_avg_end_chan_idx = mwalib_coarse_chan_range.end * fine_chans_per_coarse;
        let pre_avg_selected_frequencies = &context.metafits_context.metafits_fine_chan_freqs_hz
            [pre_avg_start_chan_idx..pre_avg_end_chan_idx];
        let pre_avg_fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz as f64;

        let avg_frequencies: Vec<f64> = pre_avg_selected_frequencies
            .chunks(avg_freq)
            .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
            .collect();
        let avg_centre_chan = avg_frequencies.len() / 2;
        let avg_centre_freq_hz = avg_frequencies[avg_centre_chan];
        let avg_fine_chan_width_hz = pre_avg_fine_chan_width_hz * avg_freq as f64;

        let obs_name = context.metafits_context.obs_name.clone();
        let field_name = match obs_name.rsplit_once("_") {
            Some((field_name, _)) => field_name,
            None => obs_name.as_str(),
        };

        let num_timesteps = (mwalib_timestep_range.len() as f64 / avg_time as f64).ceil() as usize;

        Self::new(
            path,
            num_timesteps,
            mwalib_baseline_idxs.len(),
            avg_frequencies.len(),
            start_epoch,
            avg_fine_chan_width_hz,
            avg_centre_freq_hz,
            avg_centre_chan,
            phase_centre,
            Some(field_name),
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

        // Antenna position reference frame
        hdu.write_key(&mut uvfits, "FRAME", "ITRF")?;

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

        // AIPS 117 calls this TIMESYS, but Cotter calls in TIMSYS, so we do both.
        hdu.write_key(&mut uvfits, "TIMSYS", "UTC")?;
        hdu.write_key(&mut uvfits, "TIMESYS", "UTC")?;
        hdu.write_key(&mut uvfits, "ARRNAM", "MWA")?;
        hdu.write_key(&mut uvfits, "NUMORB", 0)?; // number of orbital parameters in table
        hdu.write_key(&mut uvfits, "NOPCAL", 3)?; // Nr pol calibration values / IF(N_pcal)
        hdu.write_key(&mut uvfits, "FREQID", -1)?; // Frequency setup number
        hdu.write_key(&mut uvfits, "IATUTC", 33.0)?;

        // -> EXTVER
        // ---> in AIPS117:
        // -----> on page 12, it's "Subarray number", type I
        // -----> on page 84 onwards, all examples say "Version number of table"
        // ---> in pyuvdata is 1, presumably since we're only writing a single
        //   AIPS_AN version and someone assumed it was 1 indexed
        // ---> @derwentx: I'm pretty sure the wrong description was pasted into
        // EXTVER, and it's incorrectly being used as subarray number, when it
        // should just be the version number of the table.
        hdu.write_key(&mut uvfits, "EXTVER", 1)?;

        // -> NO_IF - Number IFs (nIF)
        // ---> in AIPS117: The value of the NO IF keyword shall specify the number of spectral
        //  windows (IFs) in the data set. In the antenna file, this controls the dimension of the
        //  polarization calibration value column.
        // ---> in Cotter, this is not used.
        // ---> since we can only deal with one spectral window at the moment,
        //  this is fixed at 1, but this would change in
        //  https://github.com/MWATelescope/Birli/issues/13
        hdu.write_key(&mut uvfits, "NO_IF", 1)?;

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
    #[allow(clippy::too_many_arguments)]
    #[inline]
    pub fn write_vis_row(
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

    /// Write visibilty and weight rows into the uvfits file from the provided
    /// [`mwalib::CorrelatorContext`].
    ///
    /// # Details
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
    /// `jones_array` a [`ndarray::Array3`] of [`Jones`] visibilities with dimensions
    /// [timestep][channel][baselines]
    ///
    /// `flag_array` a [`ndarray::Array3`] of boolean flags with dimensions
    /// identical dimensions to `jones_array`
    ///
    /// `mwalib_timestep_range` the range of timestep indices (according to mwalib)
    /// which are used in the visibility and flag arrays
    ///
    /// `mwalib_coarse_chan_range` the range of coarse channel indices (according to mwalib)
    /// which are used in the visibility and flag arrays
    ///
    /// `mwalib_baseline_idxs` the baseline indices (according to mwalib) used
    /// in the visibility and flag arrays
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
    pub fn write_jones_flags(
        &mut self,
        context: &CorrelatorContext,
        jones_array: &Array3<Jones<f32>>,
        flag_array: &Array3<bool>,
        mwalib_timestep_range: &Range<usize>,
        mwalib_coarse_chan_range: &Range<usize>,
        mwalib_baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
        draw_progress: bool,
    ) -> Result<(), IOError> {
        let num_pols = context.metafits_context.num_visibility_pols;
        let expanded_flag_array = expand_flag_array(flag_array.view(), num_pols);
        let weight_factor = get_weight_factor(context);
        let weight_array = flag_to_weight_array(expanded_flag_array.view(), weight_factor);
        self.write_vis_mwalib(
            jones_array.view(),
            weight_array.view(),
            expanded_flag_array.view(),
            context,
            mwalib_timestep_range,
            mwalib_coarse_chan_range,
            mwalib_baseline_idxs,
            avg_time,
            avg_freq,
            draw_progress,
        )
    }
}

// In order to get rid of context, would need to :
// - replace timestep range with slice of epochs
// - replace coarse chan range with...
// - replace baseline indices with ...
// - pass in tiles xyz geodetic
impl WriteableVis for UvfitsWriter {
    fn write_vis_mwalib(
        &mut self,
        jones_array: ArrayView3<Jones<f32>>,
        weight_array: ArrayView4<f32>,
        flag_array: ArrayView4<bool>,
        context: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
        avg_time: usize,
        avg_freq: usize,
        draw_progress: bool,
    ) -> Result<(), IOError> {
        let mut uvfits = self.open()?;

        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let num_sel_timesteps = timestep_range.len();
        let num_avg_timesteps = (num_sel_timesteps as f64 / avg_time as f64).ceil() as usize;

        let num_sel_coarse_chans = coarse_chan_range.len();
        let num_sel_chans = fine_chans_per_coarse * num_sel_coarse_chans;
        let num_avg_chans = (num_sel_chans as f64 / avg_freq as f64).ceil() as usize;

        let num_baselines = baseline_idxs.len();

        let num_pols = context.metafits_context.num_visibility_pols;

        let expected_dims = (num_sel_timesteps, num_sel_chans, num_baselines);
        let jones_dims = jones_array.dim();
        if jones_dims != expected_dims {
            return Err(IOError::BadArrayShape {
                argument: "jones_array".into(),
                function: "UvfitsWriter::write_vis_mwalib".into(),
                expected: format!("{:?}", expected_dims),
                received: format!("{:?}", jones_dims),
            });
        }
        let weight_dims = weight_array.dim();
        if weight_dims != (jones_dims.0, jones_dims.1, jones_dims.2, num_pols) {
            return Err(IOError::BadArrayShape {
                argument: "weight_array".into(),
                function: "UvfitsWriter::write_vis_mwalib".into(),
                expected: format!("{:?}", (jones_dims.0, jones_dims.1, jones_dims.2, 4)),
                received: format!("{:?}", weight_dims),
            });
        }

        assert_eq!(self.total_num_rows, num_avg_timesteps * num_baselines);

        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);

        // the gps start times [seconds] of all the selected timesteps
        let sel_gps_times_ms: Vec<u64> = context.timesteps[timestep_range.clone()]
            .iter()
            .map(|t| t.gps_time_ms)
            .collect();

        let sel_baselines = baseline_idxs
            .iter()
            .map(|&idx| context.metafits_context.baselines[idx].clone())
            .collect::<Vec<_>>();

        // integration time [seconds] before averaging
        let int_time_ms = context.metafits_context.corr_int_time_ms as u64;
        let avg_int_time_ms = avg_time as u64 * int_time_ms;
        let half_avg_int_time_ms = avg_int_time_ms / 2;

        let draw_target = if draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        };

        // Create a progress bar to show the writing status
        let write_progress =
            indicatif::ProgressBar::with_draw_target(self.total_num_rows as u64, draw_target);
        write_progress.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .progress_chars("=> "),
        );
        write_progress.set_message("write uv vis");

        // temporary vector to hold visiblity data, avoids heap allocation
        let mut tmp_vis: Vec<f32> = vec![0.0; 12 * num_avg_chans];

        // the result of averging a chunk
        let mut avg_jones: Jones<f32> = Jones::default();
        let mut avg_weight = Array::from_elem(4, 0.0);
        let mut avg_flag = Array::from_elem(4, false);

        for (
            gps_times_chunk_ms,
            jones_timestep_chunk,
            weight_timestep_chunk,
            flag_timestep_chunk,
        ) in izip!(
            sel_gps_times_ms.chunks(avg_time),
            jones_array.axis_chunks_iter(Axis(0), avg_time),
            weight_array.axis_chunks_iter(Axis(0), avg_time),
            flag_array.axis_chunks_iter(Axis(0), avg_time),
        ) {
            let midpoint_epoch = gps_millis_to_epoch(gps_times_chunk_ms[0] + half_avg_int_time_ms);

            let prec_info = precess_time(
                self.phase_centre,
                midpoint_epoch,
                self.array_pos.longitude_rad,
                self.array_pos.latitude_rad,
            );

            let tiles_xyz_precessed = prec_info.precess_xyz_parallel(&tiles_xyz_geod);

            for (baseline, jones_baseline_chunk, weight_baseline_chunk, flag_baseline_chunk) in izip!(
                sel_baselines.iter(),
                jones_timestep_chunk.axis_iter(Axis(2)),
                weight_timestep_chunk.axis_iter(Axis(2)),
                flag_timestep_chunk.axis_iter(Axis(2)),
            ) {
                let ant1_idx = baseline.ant1_index;
                let ant2_idx = baseline.ant2_index;

                let baseline_xyz_precessed =
                    tiles_xyz_precessed[ant1_idx] - tiles_xyz_precessed[ant2_idx];
                let uvw = UVW::from_xyz(baseline_xyz_precessed, prec_info.hadec_j2000);

                // MWA/CASA/AOFlagger visibility order is XX,XY,YX,YY
                // UVFits visibility order is XX,YY,XY,YX

                for (jones_chunk, weight_chunk, flag_chunk, vis_chunk) in izip!(
                    jones_baseline_chunk.axis_chunks_iter(Axis(1), avg_freq),
                    weight_baseline_chunk.axis_chunks_iter(Axis(1), avg_freq),
                    flag_baseline_chunk.axis_chunks_iter(Axis(1), avg_freq),
                    tmp_vis.chunks_exact_mut(12)
                ) {
                    if avg_time == 1 && avg_freq == 1 {
                        avg_jones = jones_chunk[[0, 0]];
                        avg_weight.assign(&weight_chunk.slice(s![0, 0, ..]));
                        avg_flag.assign(&flag_chunk.slice(s![0, 0, ..]));
                    } else {
                        average_chunk_for_pols_f64!(
                            jones_chunk,
                            weight_chunk,
                            flag_chunk,
                            avg_jones,
                            avg_weight,
                            avg_flag
                        );
                    }
                    izip!(avg_weight.iter_mut(), avg_flag.iter())
                        .for_each(|(w, f)| *w = if *f { -*w } else { *w });

                    // TODO: something a little prettier
                    vis_chunk[0] = avg_jones[0].re;
                    vis_chunk[1] = avg_jones[0].im;
                    vis_chunk[2] = avg_weight[0];
                    vis_chunk[3] = avg_jones[3].re;
                    vis_chunk[4] = avg_jones[3].im;
                    vis_chunk[5] = avg_weight[3];
                    vis_chunk[6] = avg_jones[1].re;
                    vis_chunk[7] = avg_jones[1].im;
                    vis_chunk[8] = avg_weight[1];
                    vis_chunk[9] = avg_jones[2].re;
                    vis_chunk[10] = avg_jones[2].im;
                    vis_chunk[11] = avg_weight[2];
                    // dbg!(&write_progress.position(), &jones_chunk, &jones_chunk.len(), &weight_chunk, &avg_jones, &avg_weight, &vis_chunk);
                }

                self.write_vis_row(
                    &mut uvfits,
                    uvw,
                    ant1_idx,
                    ant2_idx,
                    midpoint_epoch,
                    &tmp_vis,
                )?;

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
    use tempfile::NamedTempFile;

    use float_cmp::{approx_eq, F32Margin, F64Margin};

    use fitsio::{
        errors::check_status as fits_check_status,
        hdu::{FitsHdu, HduInfo},
    };

    use crate::{
        approx::assert_abs_diff_eq,
        context_to_jones_array, get_flaggable_timesteps, init_flag_array,
        marlu::{
            constants::{
                COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
            },
            mwalib::{
                _get_fits_col, _get_required_fits_key, _open_fits, _open_hdu, fits_open,
                fits_open_hdu, get_fits_col, get_required_fits_key,
            },
            time,
        },
    };

    // TODO: dedup this from lib.rs
    pub fn get_mwa_ord_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    #[allow(dead_code)]
    pub(crate) fn get_1254670392_avg_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1254670392_avg/1254670392.fixed.metafits";
        let gpufits_paths = [
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox01_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox02_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox03_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox04_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox05_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox06_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox07_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox08_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox09_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox10_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox11_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox12_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox13_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox14_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox15_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox16_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox17_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox18_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox19_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox20_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox21_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox22_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox23_00.fits",
            "tests/data/1254670392_avg/1254670392_20191009153257_gpubox24_00.fits",
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
                        assert_eq!(left_val, right_val, "mismatch for short string key {}", key,);
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
    pub(crate) use assert_short_string_keys_eq;

    macro_rules! assert_f32_string_keys_eq {
        ($keys:expr, $left_fptr:expr, $left_hdu:expr, $right_fptr:expr, $right_hdu:expr) => {
            for key in $keys {
                match (
                    get_required_fits_key!($left_fptr, &$left_hdu, key),
                    get_required_fits_key!($right_fptr, &$right_hdu, key),
                ) {
                    (Ok(left_val), Ok(right_val)) => {
                        assert!(
                            approx_eq!(f32, left_val, right_val, F32Margin::default()),
                            "mismatch for short f32 key {}. {} != {}",
                            key,
                            left_val,
                            right_val
                        );
                    }
                    (Err(err), Ok(right_val)) => {
                        panic!(
                            "unable to get left short f32 key {}. Right val={}. err={}",
                            key, right_val, err
                        );
                    }
                    (.., Err(err)) => {
                        panic!("unable to get right short f32 key {}. {}", key, err);
                    }
                }
            }
        };
    }
    pub(crate) use assert_f32_string_keys_eq;

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
            }
        };
    }
    pub(crate) use assert_table_column_descriptions_match;

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
    pub(crate) use assert_table_string_column_values_match;

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
    pub(crate) use assert_table_f64_column_values_match;

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
    pub(crate) use assert_table_vector_f64_column_values_match;

    #[allow(dead_code)]
    pub(crate) fn assert_uvfits_primary_header_eq(
        left_fptr: &mut FitsFile,
        right_fptr: &mut FitsFile,
    ) {
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

        assert_f32_string_keys_eq!(
            vec![
                "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "NAXIS4", "NAXIS5", "NAXIS6",
                "BSCALE", "PSCAL1", "PZERO1", "PSCAL2", "PZERO2", "PSCAL3", "PZERO3", "PSCAL4",
                "PZERO4", "PSCAL5", "PZERO5", "CRVAL2", "CRPIX2", "CDELT2", "CRVAL3", "CDELT3",
                "CRPIX3",
                // This is actually incorrect in Cotter, see https://github.com/MWATelescope/Birli/issues/6
                // "CRVAL4",
                "CDELT4", "CRPIX4", "CRVAL5", "CRPIX5", "CDELT5", "CRVAL6", "CRPIX6", "CDELT6",
                "EPOCH", "OBSRA",
                "OBSDEC",
                // "FIBRFACT"
                // the velocity factor of electic fields in RG-6 like coax.
                // this is incorrect in Cotter. It writes 2.0 instead of 1.204
            ],
            left_fptr,
            _left_primary_hdu,
            right_fptr,
            _right_primary_hdu
        );
    }

    #[allow(dead_code)]
    pub(crate) fn assert_uvfits_ant_header_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
        let mut _left_ant_hdu = fits_open_hdu!(left_fptr, 1).unwrap();
        let mut _right_ant_hdu = fits_open_hdu!(right_fptr, 1).unwrap();

        assert_short_string_keys_eq!(
            vec![
                "XTENSION", "TTYPE1", "TFORM1", "TTYPE2", "TFORM2", "TUNIT2", "TTYPE3", "TFORM3",
                "TTYPE4", "TFORM4", "TTYPE5", "TFORM5", "TUNIT5", "TTYPE6", "TFORM6", "TTYPE7",
                "TFORM7", "TUNIT7", "TTYPE8", "TFORM8", "TTYPE9", "TFORM9", "TTYPE10", "TFORM10",
                "TUNIT10", "TTYPE11", "TFORM11", "EXTNAME", "TIMSYS", "ARRNAM", "RDATE"
            ],
            left_fptr,
            _left_ant_hdu,
            right_fptr,
            _right_ant_hdu
        );

        assert_f32_string_keys_eq!(
            vec![
                "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS", //
                "ARRAYX", "ARRAYY", "ARRAYZ",
                // "FREQ", // This is incorrect in Cotter. See: https://github.com/MWATelescope/Birli/issues/6
                // "GSTIA0", // TODO: this is off from Cotter by about 5e-2
                "DEGPDY", "POLARX", "POLARY", "UT1UTC", "DATUTC", "NUMORB", "NOPCAL", "FREQID",
                "IATUTC",
            ],
            left_fptr,
            _left_ant_hdu,
            right_fptr,
            _right_ant_hdu
        );
    }

    #[allow(dead_code)]
    pub(crate) fn get_group_column_description(
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

    #[allow(dead_code)]
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
            }
        };
    }
    pub(crate) use assert_group_column_descriptions_match;

    #[allow(dead_code)]
    pub(crate) fn assert_uvfits_vis_table_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
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
                // ffgpve = fits_read_sel_flt
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

                // ffgpve = fits_read_sel_flt
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

    #[allow(dead_code)]
    pub(crate) fn assert_uvfits_ant_table_eq(left_fptr: &mut FitsFile, right_fptr: &mut FitsFile) {
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
    pub(crate) fn uvfits_from_mwalib_matches_cotter_header() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();
        // let tmp_uvfits_file = Path::new("tests/data/test_header.uvfits");

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = context.common_coarse_chan_indices.clone();
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            1,
            1,
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

                u.write_vis_row(
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

    /// This test ensures center frequencies are calculated correctly.
    /// See: https://github.com/MWATelescope/Birli/issues/6
    #[test]
    fn center_frequencies() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            1,
            1,
        )
        .unwrap();

        // Prepare our flagmasks with known bad antennae
        let flag_array = init_flag_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            None,
            None,
        );

        let (jones_array, flag_array) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            Some(flag_array),
            false,
        )
        .unwrap();

        u.write_jones_flags(
            &context,
            &jones_array,
            &flag_array,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            1,
            1,
            false,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();

        let expected_center_freq = 229760000.;
        let expected_fine_chan_width = 640000.;

        let birli_vis_hdu = fits_open_hdu!(&mut birli_fptr, 0).unwrap();
        let birli_vis_freq: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "CRVAL4").unwrap();
        assert_abs_diff_eq!(birli_vis_freq, expected_center_freq);
        let birli_vis_width: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "CDELT4").unwrap();
        assert_abs_diff_eq!(birli_vis_width, expected_fine_chan_width);
        let birli_ant_hdu = fits_open_hdu!(&mut birli_fptr, 1).unwrap();
        let birli_ant_freq: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "FREQ").unwrap();
        assert_abs_diff_eq!(birli_ant_freq, expected_center_freq);
    }

    /// This test ensures center frequencies are calculated correctly with frequency averaging.
    /// See: https://github.com/MWATelescope/Birli/issues/6
    #[test]
    fn avg_center_frequencies() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        assert_eq!(sel_coarse_chan_range.len(), 2);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let (avg_time, avg_freq) = (1, 2);

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            avg_time,
            avg_freq,
        )
        .unwrap();

        // Prepare our flagmasks with known bad antennae
        let flag_array = init_flag_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            None,
            None,
        );

        let (jones_array, flag_array) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            Some(flag_array),
            false,
        )
        .unwrap();

        let num_pols = context.metafits_context.num_visibility_pols;
        let flag_array = expand_flag_array(flag_array.view(), num_pols);
        let weight_factor = get_weight_factor(&context);
        let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        u.write_vis_mwalib(
            jones_array.view(),
            weight_array.view(),
            flag_array.view(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            avg_time,
            avg_freq,
            false,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();

        let expected_center_freq = (229760000. + 230400000.) / 2.;
        let expected_fine_chan_width = 1280000.;

        let birli_vis_hdu = fits_open_hdu!(&mut birli_fptr, 0).unwrap();
        let birli_vis_freq: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "CRVAL4").unwrap();
        assert_abs_diff_eq!(birli_vis_freq, expected_center_freq);
        let birli_vis_width: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "CDELT4").unwrap();
        assert_abs_diff_eq!(birli_vis_width, expected_fine_chan_width);
        let birli_ant_hdu = fits_open_hdu!(&mut birli_fptr, 1).unwrap();
        let birli_ant_freq: f64 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "FREQ").unwrap();
        assert_abs_diff_eq!(birli_ant_freq, expected_center_freq);
    }

    /// Tests for AIPS 117 compliance.
    /// See https://github.com/MWATelescope/Birli/issues/9
    /// and ftp://ftp.aoc.nrao.edu/pub/software/aips/TEXT/PUBL/AIPSMEM117.PS
    #[test]
    fn aips_117() {
        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            1,
            1,
        )
        .unwrap();

        let flag_array = init_flag_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            None,
            None,
        );

        let (jones_array, flag_array) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            Some(flag_array),
            false,
        )
        .unwrap();

        u.write_jones_flags(
            &context,
            &jones_array,
            &flag_array,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            1,
            1,
            false,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();

        // /////////// //
        // PRIMARY HDU //
        // /////////// //

        let birli_vis_hdu = fits_open_hdu!(&mut birli_fptr, 0).unwrap();

        // -> OBJECT
        let birli_vis_object: String =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "OBJECT").unwrap();
        assert_eq!(birli_vis_object, "ForA");
        // -> TELESCOP
        let birli_vis_telescop: String =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "TELESCOP").unwrap();
        assert_eq!(birli_vis_telescop, "MWA");
        // -> INSTRUME
        let birli_vis_instrume: String =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "INSTRUME").unwrap();
        assert_eq!(birli_vis_instrume, "MWA");
        // -> DATE-OBS
        // ---> in AIPS117:
        // -----> on page 12: "The value of the RDATE parameter will be the date for which the time
        //  system parameters GSTIA0, DECPDY, and IATUTC apply.
        //  If the table contains orbital parameters for orbiting antennæ, this keyword also
        //  designates the epoch for the orbital parameters.
        //  (This is copy-pasted twice)
        // -----> on page 85 onwards, all examples show YYYY-MM-DD format
        // ---> in Cotter, it is given in ISO8601 (YYYY-MM-DDTHH:mm:ss) with time fixed to 00:00:00.
        // TODO: determine whether this field should have the time.
        // let birli_vis_date_obs: String =
        // get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "DATE-OBS").unwrap();
        // assert_eq!(birli_vis_date_obs, "2017-12-01T14:54:38");

        // -> DATE-MAP - File processing date
        // ---> not in Cotter, not mandatory, so not written
        // -> BSCALE
        let birli_ant_bscale: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "BSCALE").unwrap();
        assert_abs_diff_eq!(birli_ant_bscale, 1.);
        // -> BUNIT - units,
        // ---> not in Cotter, not mandatory, so not written
        // let birli_ant_bunit: String =
        // get_required_fits_key!(&mut birli_fptr, &birli_vis_hdu, "BUNIT").unwrap();
        // -> EQUINOX - Equinox of source coordinates and uvw
        // ---> not in Cotter, not mandatory, so not written
        // -> ALTRPIX - Reference pixel for velocity
        // ---> not in Cotter, not mandatory, so not written

        // /////////// //
        // ANTENNA HDU //
        // /////////// //

        let birli_ant_hdu = fits_open_hdu!(&mut birli_fptr, 1).unwrap();

        // -> EXTVER
        // ---> in AIPS117:
        // -----> on page 12, it's "Subarray number", type I
        // -----> on page 84 onwards, all examples say "Version number of table"
        // ---> in pyuvdata is 1, presumably since we're only writing a single
        //   AIPS_AN version and someone assumed it was 1 indexed
        // ---> @derwentx: I'm pretty sure the wrong description was pasted into
        // EXTVER, and it's incorrectly being used as subarray number, when it
        // should just be the version number of the table.
        let birli_ant_extver: i32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "EXTVER").unwrap();
        assert_eq!(birli_ant_extver, 1);
        // -> ARRAY{X|Y|Z} = coordinate of array center (meters)
        let birli_ant_arrayx: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "ARRAYX").unwrap();
        assert_abs_diff_eq!(birli_ant_arrayx, -2_559_453.3, epsilon = f32::EPSILON);
        let birli_ant_arrayy: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "ARRAYY").unwrap();
        assert_abs_diff_eq!(birli_ant_arrayy, 5_095_371.5, epsilon = f32::EPSILON);
        let birli_ant_arrayz: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "ARRAYZ").unwrap();
        assert_abs_diff_eq!(birli_ant_arrayz, -2_849_056.8, epsilon = f32::EPSILON);
        // -> GSTIAO = GST at 0h on reference date (degrees)
        // ---> our value is out from cotter's by about 5e-2. Not sure if this is an issue.
        // ---> notes from @mkolopanis: this one depends on your RDATE but I'm not 100% sure what
        //   sets the RDATE for our kind of data. We just use 0 to spoof this because we don't
        //   directly use it in any of our calculations.
        let birli_ant_gstia0: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "GSTIA0").unwrap();
        assert_abs_diff_eq!(birli_ant_gstia0, 70.044_16, epsilon = 5e-2);
        // -> DEGPDY = Earth’s rotation rate (degrees/day)
        let birli_ant_degpdy: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "DEGPDY").unwrap();
        assert_abs_diff_eq!(birli_ant_degpdy, 360.985, epsilon = f32::EPSILON);
        // -> FREQ = Reference frequency (Hz)
        // ---> Cotter-calculated value was slightly incorrect because of
        //   https://github.com/MWATelescope/Birli/issues/6
        let birli_ant_freq: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "FREQ").unwrap();
        assert_abs_diff_eq!(birli_ant_freq, 229760000., epsilon = f32::EPSILON);
        // -> RDATE = Reference date
        // ---> in AIPS117:
        // -----> on page 12: "The value of the RDATE parameter will be the date for which the time
        //  system parameters GSTIA0, DECPDY, and IATUTC apply.
        //  If the table contains orbital parameters for orbiting antennæ, this keyword also
        //  designates the epoch for the orbital parameters.
        //  (This is copy-pasted twice)
        // -----> on page 85 onwards, all examples show YYYY-MM-DD format
        // ---> in Cotter, it is given in ISO8601 (YYYY-MM-DDTHH:mm:ss) with time fixed to 00:00:00.
        let birli_ant_rdate: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "RDATE").unwrap();
        assert_eq!(birli_ant_rdate, "2017-12-01T00:00:00.0");
        // -> POLAR{X|Y} = coordinate of North Pole (arc seconds)
        // ---> notes from @mkolopanis: 0 makes sense here I think. Unless there's an offset from
        //   North Pole in images.
        let birli_ant_polarx: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "POLARX").unwrap();
        assert_abs_diff_eq!(birli_ant_polarx, 0., epsilon = f32::EPSILON);
        let birli_ant_polary: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "POLARX").unwrap();
        assert_abs_diff_eq!(birli_ant_polary, 0., epsilon = f32::EPSILON);
        // -> UT1UTC = UT1 - UTC (sec)
        // ---> notes from @mkolopanis: we also use 0 here.
        let birli_ant_ut1utc: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "UT1UTC").unwrap();
        assert_abs_diff_eq!(birli_ant_ut1utc, 0., epsilon = f32::EPSILON);
        // -> DATUTC = "time system - UTC (sec)" (huh)
        // ---> notes from @mkolopanis: would 0 if your TIMESYS is UTC. We assume UTC ourselves so
        //   it is always 0
        let birli_ant_datutc: f32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "DATUTC").unwrap();
        assert_abs_diff_eq!(birli_ant_datutc, 0., epsilon = f32::EPSILON);

        // -> TIMSYS/TIMESYS "Time system"
        // ---> in AIPS117:
        // -----> on page 12 it's listed as `TIMESYS`, type A in Mandatory keywords
        // -----> on page 13 it's Time system. The TIMSYS keyword shall specify the time system used
        //        for the array. It shall either have the value ’IAT’, denoting international atomic
        //        time, or the value ’UTC’, denoting coordinated universal time. This indicates
        //        whether the zero hour for the TIME parameter in the UV DATA table is midnight IAT
        //        or midnight UTC.
        // -----> on page 87, it's `TIMSYS` in an example
        // ---> in CASA it's `TIMSYS`
        // ---> in pyuvdata: they look for both
        // ---> in Cotter it's `TIMSYS`
        // So ¯\_(ツ)_/¯
        let birli_ant_timesys: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "TIMESYS").unwrap();
        assert_eq!(birli_ant_timesys, "UTC");
        let birli_ant_timsys: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "TIMSYS").unwrap();
        assert_eq!(birli_ant_timsys, "UTC");
        // -> ARRNAM
        let birli_ant_timsys: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "ARRNAM").unwrap();
        assert_eq!(birli_ant_timsys, "MWA");
        // -> XYZHAND
        // ---> Birli assumes the station coordinates are "right handed".
        let birli_ant_frame: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "XYZHAND").unwrap();
        assert_eq!(birli_ant_frame, "RIGHT");
        // -> FRAME - Coordinate frame
        // ---> in AIPS117: The value of the FRAME keyword shall be a string that identifies the
        //  coordinate system used for antenna coordinates. At present, only one value of the FRAME
        //  keyword has been defined (’ITRF’), although ’????’ is widely used to reflect ignorance.
        // ---> notes from @mkolopanis: the "FRAME" keyword in particular in the "AIPS AN" table is
        //   important to know which frame the antenna positions are recorded in.
        let birli_ant_frame: String =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "FRAME").unwrap();
        assert_eq!(birli_ant_frame, "ITRF");
        // -> NUMORB - Number orbital parameters in table (norb)
        let birli_ant_numorb: i32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "NUMORB").unwrap();
        assert_eq!(birli_ant_numorb, 0);
        // -> NO_IF - Number IFs (nIF)
        // ---> in AIPS117: The value of the NO IF keyword shall specify the number of spectral
        //  windows (IFs) in the data set. In the antenna file, this controls the dimension of the
        //  polarization calibration value column.
        // ---> in Cotter, this is not used.
        // ---> notes from @mkolopanis: we handle NO_IF as the number of independent spectral
        //  windows. I can see how you could interpret each coarse band as its own spectral window.
        // We just view the full band as a single window since it is contiguous with all coarse
        // bands present. The idea of multiple spectral windows I think is more common for things
        // like the VLA when you observe in multiple bands that do not form a contiguous frequency
        // band if concatenated.
        let birli_ant_no_if: i32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "NO_IF").unwrap();
        assert_eq!(birli_ant_no_if, 1);
        // -> NOPCAL - Number of polarization calibration values / IF (npcal)
        let birli_ant_nopcal: i32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "NOPCAL").unwrap();
        assert_eq!(birli_ant_nopcal, 3);
        // -> POLTYPE - Type of polarization calibration
        // ---> If the table contains information about the polarization characteristics of
        //  the feeds, then the feed parametrization that is used shall be indicated by the value of
        //  the POLTYPE keyword, as given in Table 9.
        // ---> not used.
        // let birli_ant_poltype: i32 =
        //     get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "POLTYPE").unwrap();

        // -> FREQID - Frequency setup number
        let birli_ant_freqid: i32 =
            get_required_fits_key!(&mut birli_fptr, &birli_ant_hdu, "FREQID").unwrap();
        assert_eq!(birli_ant_freqid, -1);
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use super::*;
    use tempfile::NamedTempFile;
    use tests::{
        assert_uvfits_ant_header_eq, assert_uvfits_ant_table_eq, assert_uvfits_primary_header_eq,
        assert_uvfits_vis_table_eq, get_1254670392_avg_context, get_mwa_ord_context,
    };

    use crate::{
        context_to_jones_array,
        flags::flag_jones_array_existing,
        get_antenna_flags, get_baseline_flags, get_flaggable_timesteps, init_flag_array,
        marlu::{
            constants::{
                COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
            },
            mwalib::{_open_fits, fits_open},
        },
    };
    use aoflagger_sys::cxx_aoflagger_new;

    #[test]
    fn uvfits_tables_from_mwalib_matches_cotter() {
        use crate::flags::get_baseline_flags;

        let context = get_mwa_ord_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            1,
            1,
        )
        .unwrap();

        let aoflagger = unsafe { cxx_aoflagger_new() };

        // Prepare our flagmasks with known bad antennae
        let flag_array = init_flag_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            None,
            Some(&get_baseline_flags(&context, &get_antenna_flags(&context))),
        );

        let (jones_array, flag_array) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            Some(flag_array),
            false,
        )
        .unwrap();

        let strategy_path = &aoflagger.FindStrategyFileMWA();

        let flag_array = flag_jones_array_existing(
            &aoflagger,
            strategy_path,
            &jones_array,
            Some(flag_array),
            true,
            false,
        );

        u.write_jones_flags(
            &context,
            &jones_array,
            &flag_array,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            1,
            1,
            false,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        let cotter_uvfits_path = Path::new("tests/data/1196175296_mwa_ord/1196175296.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_primary_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_vis_table_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_table_eq(&mut birli_fptr, &mut cotter_fptr);
    }

    #[test]
    fn uvfits_tables_from_mwalib_matches_cotter_avg_4s_160khz() {
        let context = get_1254670392_avg_context();

        let tmp_uvfits_file = NamedTempFile::new().unwrap();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
        let sel_baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();

        let array_pos = Some(LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        });

        let (avg_time, avg_freq) = (2, 4);

        let mut u = UvfitsWriter::from_mwalib(
            tmp_uvfits_file.path(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            array_pos,
            None,
            avg_time,
            avg_freq,
        )
        .unwrap();

        let aoflagger = unsafe { cxx_aoflagger_new() };

        // Prepare our flagmasks with known bad antennae
        let flag_array = init_flag_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            None,
            Some(&get_baseline_flags(&context, &get_antenna_flags(&context))),
        );

        let (jones_array, flag_array) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            Some(flag_array),
            false,
        )
        .unwrap();

        let strategy_path = &aoflagger.FindStrategyFileMWA();

        let flag_array = flag_jones_array_existing(
            &aoflagger,
            strategy_path,
            &jones_array,
            Some(flag_array),
            true,
            false,
        );

        let num_pols = context.metafits_context.num_visibility_pols;
        let flag_array = expand_flag_array(flag_array.view(), num_pols);
        let weight_factor = get_weight_factor(&context);
        let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        u.write_vis_mwalib(
            jones_array.view(),
            weight_array.view(),
            flag_array.view(),
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            &sel_baseline_idxs,
            avg_time,
            avg_freq,
            false,
        )
        .unwrap();

        u.write_ants_from_mwalib(&context.metafits_context).unwrap();

        let cotter_uvfits_path =
            Path::new("tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.uvfits");

        let mut birli_fptr = fits_open!(&tmp_uvfits_file.path()).unwrap();
        let mut cotter_fptr = fits_open!(&cotter_uvfits_path).unwrap();

        assert_uvfits_primary_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_vis_table_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_header_eq(&mut birli_fptr, &mut cotter_fptr);
        assert_uvfits_ant_table_eq(&mut birli_fptr, &mut cotter_fptr);
    }
}
