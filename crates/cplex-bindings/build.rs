use std::{collections::BTreeMap, env, path::PathBuf};

use bindgen::callbacks::{MacroParsingBehavior, ParseCallbacks};

const FP_EXCLUDE_MACROS: &[&str] = &[
    "FP_NAN",
    "FP_INFINITE",
    "FP_ZERO",
    "FP_SUBNORMAL",
    "FP_NORMAL",
];

// TODO: Make these constants configurable.
// Default CPLEX installation path
const CPLEX_PATH: &str = "/opt/ibm/ILOG";
// The version I have installed
const CPLEX_VERSION: &str = "CPLEX_Studio221";

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
        "cargo:rustc-link-search={CPLEX_PATH}/{CPLEX_VERSION}/cplex/lib/x86-64_linux/static_pic/"
    );

    // Links with with /opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/libcplex.a
    println!("cargo:rustc-link-lib=static=cplex");

    let bindings = bindgen::Builder::default()
        .header("cplex_wrapper.h")
        .clang_arg(format!(
            "-I{CPLEX_PATH}/{CPLEX_VERSION}/cplex/include/ilcplex/"
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
