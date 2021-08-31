extern crate cxx_build;
use std::env;

fn main() {
    println!("cargo:rustc-cfg=verbose");

    if cfg!(feature = "aoflagger") {
        match env::var("DOCS_RS").as_deref() {
            Ok("1") => (),
            _ => {
                cxx_build::bridge("src/cxx_aoflagger.rs")
                    .flag_if_supported("-std=c++11")
                    .flag_if_supported("-Wno-nonportable-include-path")
                    .flag_if_supported("-Wno-unused-parameter")
                    .file("src/cxx_aoflagger.cc")
                    .compile("cxx_aoflagger");
                // Tell cargo to tell rustc to link the aoflagger shared library.
                println!("cargo:rustc-link-lib=dylib=aoflagger");
            }
        }
        println!("cargo:rerun-if-changed=src/cxx_aoflagger.rs");
        println!("cargo:rerun-if-changed=src/cxx_aoflagger.cc");
    }

    println!("cargo:rerun-if-changed=build.rs");
}
