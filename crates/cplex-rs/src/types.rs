//! Various miscellaneous types used by the CPLEX library.

use cplex_bindings as cplex;

pub type CplexInt = std::os::raw::c_int;
pub type CplexStatus = std::os::raw::c_int;
pub type CplexDouble = std::os::raw::c_double;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Xctype {
    Continuous = cplex::CPX_CONTINUOUS,
    Binary = cplex::CPX_BINARY,
    Integer = cplex::CPX_INTEGER,
    Semicount = cplex::CPX_SEMICONT,
    Semiint = cplex::CPX_SEMIINT,
}

// A column in a cplex problem (ie a variable)
// TODO: What's the xctype?
pub struct CplexCol {
    pub obj: CplexDouble,
    pub lb: CplexDouble,
    pub ub: CplexDouble,
    pub xctype: Xctype,
    // Also, what's the semantics of this in the C code?
    // Ie does it copy from my buffer or does it assume that the string is static?
    // Apparently cplex wants *mut *mut c_char, so I guess CString is out because it's immutable.
    pub name: String,
}

// A row in a cplex problem (ie a constraint)
pub struct CplexRow {
    // The sense of the row, including the right hand sides
    pub sense: RowSense,
    // A pair of column index and coefficient
    pub coeffs: Vec<(CplexInt, CplexDouble)>,
    pub name: String,
}

#[derive(Debug)]
pub struct CplexSolution {
    pub solstat: CplexInt,
    pub objval: CplexDouble,
    pub x: Vec<CplexDouble>,
    pub pi: Vec<CplexDouble>,
    pub slack: Vec<CplexDouble>,
    pub dj: Vec<CplexDouble>,
}

#[derive(Debug, Clone, Copy)]
pub enum RowSense {
    Leq(CplexDouble),
    Eq(CplexDouble),
    Geq(CplexDouble),
    Range(CplexDouble, CplexDouble),
}

impl RowSense {
    /// Cplex wants the sense passed as a char.
    pub fn as_char(&self) -> i8 {
        match self {
            RowSense::Leq(_) => b'L' as i8,
            RowSense::Eq(_) => b'E' as i8,
            RowSense::Geq(_) => b'G' as i8,
            RowSense::Range(_, _) => b'R' as i8,
        }
    }

    pub fn rhs(&self) -> CplexDouble {
        match self {
            RowSense::Leq(rhs) => *rhs,
            RowSense::Eq(rhs) => *rhs,
            RowSense::Geq(rhs) => *rhs,
            RowSense::Range(low, _high) => *low,
        }
    }

    // Returns None if the sense is not a range.
    pub fn range(&self) -> Option<CplexDouble> {
        match self {
            RowSense::Leq(_) => None,
            RowSense::Eq(_) => None,
            RowSense::Geq(_) => None,
            RowSense::Range(low, high) => Some(high - low),
        }
    }
}
