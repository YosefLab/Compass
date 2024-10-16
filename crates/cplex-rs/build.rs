use std::{env, path::PathBuf};

pub fn main() {
    /*let mut cpx_param_regex = regex::Regex::new(r"CPXPARAM_([A-Z_]+)").unwrap();

    // Find the cplex_bindings.rs file
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    let contents = std::fs::read_to_string(out_path.join("cplex_bindings.rs")).unwrap();

    let mut cpx_params = Vec::new();
    for cap in cpx_param_regex.captures_iter(&contents) {
        cpx_params.push(cap[1].to_string());
    }

    let mut cpx_param_enum = String::new();
    cpx_param_enum.push_str("#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]\n");
    cpx_param_enum.push_str("pub enum CplexParam {\n");
    for param in cpx_params {
        cpx_param_enum.push_str(&format!("    {},\n", param));
    }
    cpx_param_enum.push_str("}\n");

    // Write the enum to a file
    std::fs::write(out_path.join("cplex_params.rs"), cpx_param_enum).unwrap();*/
}
