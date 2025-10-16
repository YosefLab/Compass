use std::{collections::BTreeMap, env, path::PathBuf};

use anyhow::Context;
use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

const FP_EXCLUDE_MACROS: &[&str] = &[
    "FP_NAN",
    "FP_INFINITE",
    "FP_ZERO",
    "FP_SUBNORMAL",
    "FP_NORMAL",
];

// Default CPLEX installation path
const DEFAULT_CPLEX_PATH: &str = "/opt/ibm/ILOG";

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

fn get_cplex_version(cplex_path: &str) -> Result<String, anyhow::Error> {
    const CPLEX_DIR_PREFIX: &str = "CPLEX_Studio";
    if let Ok(ver) = std::env::var("CPLEX_VERSION") {
        println!("Using CPLEX_VERSION from environment: {ver}");
        return Ok(ver);
    }
    println!("Searching for CPLEX version in {cplex_path}");
    let mut version: Option<String> = None;
    let entries = std::fs::read_dir(cplex_path).with_context(|| format!("Reading from {cplex_path}"))?;
    for entry in entries {
        let entry = entry.context("Failed to read directory entry")?;
        let file_os_name = entry.file_name();
        let file_name = file_os_name.to_string_lossy();
        if file_name.contains(CPLEX_DIR_PREFIX) {
            println!("Found CPLEX version {file_name}");
            match &version {
                None => version = Some(file_name.into_owned()),
                Some(prev) => {
                    // Prefer any non-community edition
                    // Prefer later versions
                    let prev_is_community = prev.contains("Community");
                    let curr_is_community = file_name.contains("Community");
                    match (prev_is_community, curr_is_community) {
                        (true, false) => version = Some(file_name.into_owned()),
                        (false, true) => {}
                        _ => {
                            // Both are community or both are non-community
                            // Prefer later version, based on lexicographical order.
                            // E.G. CPLEX_Studio2212 > CPLEX_Studio221
                            if file_name.as_ref() > prev.as_ref() {
                                version = Some(file_name.into_owned());
                            }
                        }
                    }
                }
            }
        }
    }
    match version {
        Some(version) => Ok(version),
        None => return Err(anyhow::anyhow!("Could not find CPLEX version in {cplex_path}")),
    }
}

pub fn main() {
    let cplex_path = std::env::var("CPLEX_PATH").unwrap_or_else(|_| DEFAULT_CPLEX_PATH.to_string());
    // Search for libraries in cplex path, or use the environent variable if set
    let cplex_version = get_cplex_version(&cplex_path).unwrap();
    println!("Using CPLEX version: {cplex_version}");
    // Prefer CplexStudio
    println!(
        "cargo:rustc-link-search={cplex_path}/{cplex_version}/cplex/lib/x86-64_linux/static_pic/"
    );

    // Links with with /opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/libcplex.a
    println!("cargo:rustc-link-lib=static=cplex");

    let bindings = bindgen::Builder::default()
        .header("cplex_wrapper.h")
        .clang_arg(format!(
            "-I{cplex_path}/{cplex_version}/cplex/include/ilcplex/"
        ))
        .parse_callbacks(Box::new(FpMacroExcluder))
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    let bindings_text = bindings.to_string();
    std::fs::write(out_path.join("cplex_bindings.rs"), &bindings_text).unwrap();

    // Parse the CPXPARAM_ constants into an enum
    let cpx_param_regex = regex::Regex::new(r"CPXPARAM_([A-Za-z_]+): u32 = ([0-9]+)").unwrap();

    let mut cpx_params = BTreeMap::new();
    for cap in cpx_param_regex.captures_iter(&bindings_text) {
        // Prefer later values, for some reason certain parameters are defined as
        // #define CPXPARAM_parametername 1510
        // #define CPXPARAM_ParameterName 1510
        cpx_params.insert(cap[2].parse::<u32>().unwrap(), cap[1].to_string());
    }

    let mut cpx_param_enum = String::new();
    cpx_param_enum.push_str("#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]\n");
    cpx_param_enum.push_str("#[repr(u32)]\n");
    cpx_param_enum.push_str("pub enum CplexParam {\n");
    for (_val, param) in cpx_params {
        cpx_param_enum.push_str(&format!("    {param} = CPXPARAM_{param},\n"));
    }
    cpx_param_enum.push_str("}\n");

    // Write the enum to a file
    std::fs::write(out_path.join("cplex_params.rs"), cpx_param_enum).unwrap();
}
