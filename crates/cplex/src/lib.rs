// Identifiers for the cplex library are generally not compliant with Rust's naming conventions.
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
// Also, it has extern blocks with u128. I thought u128 alignment got fixed recently though.
#![allow(improper_ctypes)]

include!(concat!(env!("OUT_DIR"), "/cplex_bindings.rs"));