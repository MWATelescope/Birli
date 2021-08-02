// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Handling of Earth Coordinates (Latitude/Longitude/Height)

use erfa_sys::{ERFA_GRS80, ERFA_WGS72, ERFA_WGS84};

use super::{error::ErfaError, XyzGeocentric};

#[derive(Clone, Debug)]
pub struct LatLng {
    /// Longitude \[radians\]
    pub longitude_rad: f64,
    /// Latitude \[radians\]
    pub latitude_rad: f64,
    /// Height above sea level \[meters\]
    pub height_metres: f64,
}

pub enum Ellipsoid {
    WGS84 = ERFA_WGS84 as isize,
    GRS80 = ERFA_GRS80 as isize,
    WGS72 = ERFA_WGS72 as isize,
}

impl LatLng {
    pub fn to_geocentric(self, ellipsoid: Ellipsoid) -> Result<XyzGeocentric, ErfaError> {
        let mut geocentric_vector: [f64; 3] = [0.0; 3];
        let status = unsafe {
            erfa_sys::eraGd2gc(
                ellipsoid as i32,               // ellipsoid identifier (Note 1)
                self.longitude_rad,             // longitude (radians, east +ve)
                self.latitude_rad,              // latitude (geodetic, radians, Note 3)
                self.height_metres,             // height above ellipsoid (geodetic, Notes 2,3)
                geocentric_vector.as_mut_ptr(), // geocentric vector (Note 2)
            )
        };
        if status != 0 {
            return Err(ErfaError {
                source_file: file!(),
                source_line: line!(),
                status,
                function: "eraGd2gc",
            });
        }
        Ok(XyzGeocentric {
            x: geocentric_vector[0],
            y: geocentric_vector[1],
            z: geocentric_vector[2],
        })
    }

    pub fn to_geocentric_wgs84(self) -> Result<XyzGeocentric, ErfaError> {
        self.to_geocentric(Ellipsoid::WGS84)
    }
}
