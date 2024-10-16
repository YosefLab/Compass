//! Example LP stolen from the CPLEX example lpex1.c
use cplex_bindings::{self as cplex, params::CplexParam};

use crate::{
    env::CplexEnv,
    problem::CplexLp,
    types::{CplexCol, CplexInt, CplexRow, RowSense, Xctype},
};

#[test]
pub fn example_linear_program() {
    let mut env = CplexEnv::new().unwrap();

    // Turn on output and data checking
    env.set_param(CplexParam::ScreenOutput, cplex::CPX_ON as CplexInt)
        .unwrap();

    env.set_param(
        CplexParam::Read_DataCheck,
        cplex::CPX_DATACHECK_WARN as CplexInt,
    )
    .unwrap();

    println!("Set parameters");

    let mut lp = CplexLp::new(&mut env, "lpex1").unwrap();

    let cols = [
        CplexCol {
            obj: 1.0,
            lb: 0.0,
            ub: 40.0,
            xctype: Xctype::Continuous,
            name: "x1".to_string(),
        },
        CplexCol {
            obj: 2.0,
            lb: 0.0,
            ub: cplex::CPX_INFBOUND_d,
            xctype: Xctype::Continuous,
            name: "x2".to_string(),
        },
        CplexCol {
            obj: 3.0,
            lb: 0.0,
            ub: cplex::CPX_INFBOUND_d,
            xctype: Xctype::Continuous,
            name: "x3".to_string(),
        },
    ];

    lp.set_objective_sense(true).unwrap();
    lp.new_cols(&cols).unwrap();

    let rows = [
        CplexRow {
            sense: RowSense::Leq(20.0),
            coeffs: vec![(0, -1.0), (1, 1.0), (2, 1.0)],
            name: "c1".to_string(),
        },
        CplexRow {
            sense: RowSense::Leq(30.0),
            coeffs: vec![(0, 1.0), (1, -3.0), (2, 1.0)],
            name: "c2".to_string(),
        },
    ];
    lp.add_rows(&rows).unwrap();

    let num_cols = lp.num_cols();
    let num_rows = lp.num_rows();
    println!("Number of columns: {num_cols}, number of rows: {num_rows}");

    lp.opt().unwrap();
    println!("Problem solved");

    let sol = lp.get_solution().unwrap();
    println!("Solution status: {sol:#?}");
}
