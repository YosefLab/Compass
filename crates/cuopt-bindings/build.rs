use std::{env, path::PathBuf};

use anyhow::{Context, anyhow};
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

const FP_EXCLUDE_MACROS: &[&str] = &[
    "FP_NAN",
    "FP_INFINITE",
    "FP_ZERO",
    "FP_SUBNORMAL",
    "FP_NORMAL",
];

// Bindgen has issues with anonymous enum and define colliding
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
    main_inner().unwrap();
}

#[derive(Debug)]
#[expect(non_camel_case_types, reason = "Matches library name")]
struct cuOpt {
    lib_dir: PathBuf,
    include_dir: PathBuf,
}

fn find_cuopt_conda() -> Result<cuOpt, anyhow::Error> {
    if let Some(conda_prefix) = env::var_os("CONDA_PREFIX") {
        println!("Found CONDA_PREFIX at {:?}", conda_prefix);
        let conda_path = std::path::PathBuf::from(conda_prefix);
        let cuopt_lib_path = conda_path.join("lib");
        let cuopt_include_path = conda_path
            .join("include");
        let lib_exists = cuopt_lib_path.exists();
        let include_exists = cuopt_include_path.exists();
        if lib_exists && include_exists {
            return Ok(cuOpt {
                lib_dir: cuopt_lib_path,
                include_dir: cuopt_include_path,
            });
        } else {
            return Err(anyhow!(
                "cuOpt conda package in unexpected format. Expected to find folders:
{} was found: {lib_exists};
{} was found: {include_exists};
This code was tested with cuOpt conda package version 25.10.00, build cuda13_251014_99e549ce. Compare your conda list cuopt",
                cuopt_lib_path.display(),
                cuopt_include_path.display(),
            ));
        }
    } else {
        return Err(anyhow!("CONDA_PREFIX environment variable not set"));
    }
}

pub fn main_inner() -> Result<(), anyhow::Error> {
    let cuopt = find_cuopt_conda().context("Failed to find CUDA")?;
    println!("Found cuOpt library at {:#?}", cuopt);

    // Expecting to find libcuopt.so
    println!("cargo:rustc-link-search={}", cuopt.lib_dir.display());
    println!("cargo:rustc-link-lib=dylib=cuopt");

    let bindings = bindgen::Builder::default()
        .header("cuopt_wrapper.h")
        .clang_arg(format!("-I{}", cuopt.include_dir.display()))
        .parse_callbacks(Box::new(FpMacroExcluder))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .context("Generating bindings")?;

    let out_path =
        PathBuf::from(env::var("OUT_DIR").context("Getting $OUT_DIR environment variable")?)
            .join("bindings.rs");
    bindings
        .write_to_file(&out_path)
        .with_context(|| format!("Writing bindings to {}", out_path.display()))?;

    Ok(())
}