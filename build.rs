// stolen from hyperdrive
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

fn main() {
    // Gather build time info
    built::write_built_file().expect("Failed to acquire build-time information");
}
