extern crate cxx_build;

fn main() {
    println!("cargo:rustc-cfg=verbose");

    cxx_build::bridge("src/cxx_aoflagger.rs")
        .flag_if_supported("-std=c++11")
        .flag_if_supported("-Wno-nonportable-include-path")
        .flag_if_supported("-Wno-unused-parameter")
        // .include("include/cxx_aoflagger.h")
        .file("src/cxx_aoflagger.cc")
        .compile("cxx_aoflagger");

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=src/cxx_aoflagger.rs");
    println!("cargo:rerun-if-changed=src/cxx_aoflagger.cc");
    println!("cargo:rerun-if-changed=include/cxx_aoflagger.h");

    // Tell cargo to tell rustc to link the aoflagger shared library.
    println!("cargo:rustc-link-lib=dylib=aoflagger");
}
