use std::{env, path::PathBuf};

use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

const FP_EXCLUDE_MACROS: &[&str] = &[
    "FP_NAN", "FP_INFINITE", "FP_ZERO", "FP_SUBNORMAL", "FP_NORMAL"
];

#[derive(Debug)]
pub struct FpMacroExcluder;

impl ParseCallbacks for FpMacroExcluder {
    fn will_parse_macro(&self, name: &str) -> MacroParsingBehavior {
        if FP_EXCLUDE_MACROS.contains(&name) {
            MacroParsingBehavior::Ignore
        } else {
            MacroParsingBehavior::Default
        }
    }
}

pub fn main() {
    println!(
        "cargo:rustc-link-search=/opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/"
    );

    // Links with with /opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/libcplex.a
    println!("cargo:rustc-link-lib=static=cplex"); 

    let bindings = bindgen::Builder::default()
        .header("cplex_wrapper.h")
        .clang_arg("-I/opt/ibm/ILOG/CPLEX_Studio221/cplex/include/ilcplex/")
        .parse_callbacks(Box::new(FpMacroExcluder))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

        let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
        bindings
            .write_to_file(out_path.join("cplex_bindings.rs"))
            .expect("Couldn't write bindings!");
}
