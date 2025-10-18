use cuopt_bindings::*;
use std::{ffi::c_char, ptr};

/// Test simple LP problem
/// Solve the following LP:
///   minimize -0.2*x1 + 0.1*x2
///   subject to:
///   3.0*x1 + 4.0*x2 <= 5.4
///   2.7*x1 + 10.1*x2 <= 4.9
///   x1, x2 >= 0
pub fn main() {
    let mut problem: cuOptOptimizationProblem = ptr::null_mut();
    let mut settings: cuOptSolverSettings = ptr::null_mut();
    let mut solution: cuOptSolution = ptr::null_mut();

    let num_variables: cuopt_int_t = 2;
    let num_constraints: cuopt_int_t = 2;

    // CSR format constraint matrix
    // https://docs.nvidia.com/nvpl/latest/sparse/storage_format/sparse_matrix.html#compressed-sparse-row-csr
    // From the constraints:
    // 3.0*x1 + 4.0*x2 <= 5.4
    // 2.7*x1 + 10.1*x2 <= 4.9
    let row_offsets: [cuopt_int_t; 3] = [0, 2, 4];
    let column_indices: [cuopt_int_t; 4] = [0, 1, 0, 1];

    let values: [cuopt_float_t; 4] = [3.0, 4.0, 2.7, 10.1];

    // Objective coefficients
    // From the objective function: minimize -0.2*x1 + 0.1*x2
    let objective_coefficients: [cuopt_float_t; 2] = [-0.2, 0.1];

    // Constraint bounds
    // From the constraints:
    // 3.0*x1 + 4.0*x2 <= 5.4
    // 2.7*x1 + 10.1*x2 <= 4.9
    let constraint_upper_bounds: [cuopt_float_t; 2] = [5.4, 4.9];
    let constraint_lower_bounds: [cuopt_float_t; 2] = [f64::NEG_INFINITY, f64::NEG_INFINITY];

    // Variable bounds
    // From the constraints: x1, x2 >= 0
    let var_lower_bounds: [cuopt_float_t; 2] = [0.0, 0.0];
    let var_upper_bounds: [cuopt_float_t; 2] = [f64::INFINITY, f64::INFINITY];

    // Variable types (continuous)
    let variable_types: [c_char; 2] = [CUOPT_CONTINUOUS as c_char, CUOPT_CONTINUOUS as c_char];

    println!("Creating and solving simple LP problem...");

    unsafe {
        // Create the problem
        let mut status = cuOptCreateRangedProblem(
            num_constraints,
            num_variables,
            CUOPT_MINIMIZE,
            0.0, // objective offset
            objective_coefficients.as_ptr(),
            row_offsets.as_ptr(),
            column_indices.as_ptr(),
            values.as_ptr(),
            constraint_lower_bounds.as_ptr(),
            constraint_upper_bounds.as_ptr(),
            var_lower_bounds.as_ptr(),
            var_upper_bounds.as_ptr(),
            variable_types.as_ptr(),
            &mut problem,
        );
        if status != CUOPT_SUCCESS {
            eprintln!("Error creating problem: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        // Create solver settings
        status = cuOptCreateSolverSettings(&mut settings);
        if status != CUOPT_SUCCESS {
            eprintln!("Error creating solver settings: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        // Set solver parameters
        status = cuOptSetFloatParameter(settings, CUOPT_ABSOLUTE_PRIMAL_TOLERANCE.as_ptr(), 0.0001);
        if status != CUOPT_SUCCESS {
            eprintln!("Error setting optimality tolerance: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        // Solve the problem
        status = cuOptSolve(problem, settings, &mut solution);
        if status != CUOPT_SUCCESS {
            eprintln!("Error solving problem: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        // Get solution information
        let mut time: cuopt_float_t = 0.0;
        status = cuOptGetSolveTime(solution, &mut time);
        if status != CUOPT_SUCCESS {
            eprintln!("Error getting solve time: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        let mut termination_status: cuopt_int_t = 0;
        status = cuOptGetTerminationStatus(solution, &mut termination_status);
        if status != CUOPT_SUCCESS {
            eprintln!("Error getting termination status: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        let mut objective_value: cuopt_float_t = 0.0;
        status = cuOptGetObjectiveValue(solution, &mut objective_value);
        if status != CUOPT_SUCCESS {
            eprintln!("Error getting objective value: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        // Print results
        println!("\nResults:");
        println!("--------");
        println!("Termination status: {}", termination_status);
        println!("Solve time: {} seconds", time);
        println!("Objective value: {}", objective_value);

        // Get and print solution variables
        let mut solution_values: Vec<cuopt_float_t> = vec![0.0; num_variables as usize];
        status = cuOptGetPrimalSolution(solution, solution_values.as_mut_ptr());
        if status != CUOPT_SUCCESS {
            eprintln!("Error getting solution values: {}", status);
            cleanup(problem, settings, solution);
            return;
        }

        println!("\nPrimal Solution: Solution variables");
        for (i, &value) in solution_values.iter().enumerate() {
            println!("x{} = {}", i + 1, value);
        }

        cleanup(problem, settings, solution);
    }
}

unsafe fn cleanup(
    mut problem: cuOptOptimizationProblem,
    mut settings: cuOptSolverSettings,
    mut solution: cuOptSolution,
) {
    unsafe {
        if !problem.is_null() {
            cuOptDestroyProblem(&mut problem);
        }
        if !settings.is_null() {
            cuOptDestroySolverSettings(&mut settings);
        }
        if !solution.is_null() {
            cuOptDestroySolution(&mut solution);
        }
    }
}
