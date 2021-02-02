use clap::{crate_authors, crate_description, crate_name, crate_version, App};
use std::env;

fn main() {
    App::new(crate_name!())
        .version(crate_version!())
        .author(crate_authors!())
        .about(crate_description!())
        .get_matches();
}

#[cfg(test)]
mod tests {
    use assert_cli;
    use std::env;

    #[test]
    fn main_with_version() {
        let pkg_name = env!("CARGO_PKG_NAME");
        let pkg_version = env!("CARGO_PKG_VERSION");
        assert_cli::Assert::main_binary()
            .with_args(&["--version"])
            .succeeds()
            .stdout()
            .is(&format!("{} {}\n", pkg_name, pkg_version)[..])
            .unwrap();
    }
}
