// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Useful constants.
//!
//! All constants *must* be double precision.
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

// cotter's constants. Useful for being more precise when converting geocentric
// XYZ to geodetic XYZ!
/// cotter's MWA longitude on Earth in radians. Use [MWA_LONG_RAD] unless you
/// know what you're doing.
pub const COTTER_MWA_LONGITUDE_RADIANS: f64 = 2.0362897754687257;
/// cotter's MWA latitude on Earth in radians. Use [MWA_LAT_RAD] unless you know
/// what you're doing.
pub const COTTER_MWA_LATITUDE_RADIANS: f64 = -0.46606083776035967;
/// cotter's MWA altitude in metres. Use [MWA_HEIGHT_M] unless you know what
/// you're doing.
pub const COTTER_MWA_HEIGHT_METRES: f64 = 377.0;

/// This is the number of seconds from 1900 Jan 1 and 1980 Jan 5. The GPS epoch
/// is 1980 Jan 5, but `hifitime` uses 1900 for everything; subtracting this
/// number from the result of `hifitime::Epoch::as_gpst_seconds` gives the
/// expected GPS time.
pub const HIFITIME_GPS_FACTOR: f64 =
    hifitime::SECONDS_PER_YEAR * 80.0 + hifitime::SECONDS_PER_DAY * 4.0;