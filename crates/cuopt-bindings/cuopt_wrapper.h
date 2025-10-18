#include <cuopt/linear_programming/cuopt_c.h>

// bindgen ends up not outputting a constant for CUOPT_INFINITY
// probably related to https://github.com/rust-lang/rust-bindgen/issues/2426
// Recommened to just use cuopt_float_t::INFINITY from the rust side.
const cuopt_float_t CUOPT_INFINITY_d = CUOPT_INFINITY;

// bindgen outputs these types 
const char CUOPT_LESS_THAN_d = CUOPT_LESS_THAN;
const char CUOPT_GREATER_THAN_d = CUOPT_GREATER_THAN;
const char CUOPT_EQUAL_d = CUOPT_EQUAL;
const char CUOPT_CONTINUOUS_d = CUOPT_CONTINUOUS;
const char CUOPT_INTEGER_d = CUOPT_INTEGER;