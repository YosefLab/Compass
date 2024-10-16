use std::{
    ffi::{CString},
    ptr::{addr_of_mut},
};

use cplex_bindings as cplex;

pub fn main() {
    println!("Hello, world!");

    unsafe {
        let mut status: std::os::raw::c_int = 0;
        let mut env = CplexEnv::new().unwrap();

        // Turn out output and data checking
        status = cplex::CPXsetintparam(
            env.inner,
            cplex::CPXPARAM_ScreenOutput as CplexInt,
            cplex::CPX_ON as CplexInt,
        );
        if status != 0 {
            let err = CplexError::new(&env, status);
            panic!("CPXopenCPLEX failed with error {err:?}");
        }

        status = cplex::CPXsetintparam(
            env.inner,
            cplex::CPXPARAM_Read_DataCheck as CplexInt,
            cplex::CPX_DATACHECK_WARN as CplexInt,
        );
        if status != 0 {
            let err = CplexError::new(&env, status);
            panic!("CPXopenCPLEX failed with error {err:?}");
        }

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

    println!("Goodbye, world!");
}

pub type CplexInt = std::os::raw::c_int;
pub type CplexStatus = std::os::raw::c_int;
pub type CplexDouble = std::os::raw::c_double;

pub struct CplexEnv {
    inner: cplex::CPXENVptr,
}

pub struct CplexLp<'a> {
    inner: cplex::CPXLPptr,
    // TODO: Arc it?
    env: &'a CplexEnv,
}

#[derive(Debug)]
pub struct CplexError {
    status: CplexStatus,
    msg: String,
}

// A column in a cplex problem (ie a variable)
// TODO: What's the xctype?
pub struct CplexCol {
    obj: CplexDouble,
    lb: CplexDouble,
    ub: CplexDouble,
    xctype: Xctype,
    // Also, what's the semantics of this in the C code?
    // Ie does it copy from my buffer or does it assume that the string is static?
    // Apparently cplex wants *mut *mut c_char, so I guess CString is out because it's immutable.
    name: String,
}

// A row in a cplex problem (ie a constraint)
pub struct CplexRow {
    // The sense of the row, including the right hand sides
    sense: RowSense,
    // A pair of column index and coefficient
    coeffs: Vec<(CplexInt, CplexDouble)>,
    name: String,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Xctype {
    Continuous = cplex::CPX_CONTINUOUS,
    Binary = cplex::CPX_BINARY,
    Integer = cplex::CPX_INTEGER,
    Semicount = cplex::CPX_SEMICONT,
    Semiint = cplex::CPX_SEMIINT,
}

impl CplexEnv {
    pub fn new() -> Result<Self, CplexError> {
        let mut status: CplexStatus = 0;
        let inner = unsafe {
            let ptr = cplex::CPXopenCPLEX(addr_of_mut!(status));
            if ptr.is_null() {
                return Err(CplexError::new_ptr(ptr, status));
            }
            ptr
        };
        Ok(Self { inner })
    }
}

impl Drop for CplexEnv {
    fn drop(&mut self) {
        println!("Dropping CplexEnv");
        unsafe {
            let status = cplex::CPXcloseCPLEX(addr_of_mut!(self.inner));
            if status != 0 {
                let err = CplexError::new_ptr(self.inner, status);
                eprintln!("CPXcloseCPLEX failed with error {err:?}");
            }
        }
    }
}

impl<'a> CplexLp<'a> {
    pub fn new(env: &'a CplexEnv, name: &str) -> Result<Self, CplexError> {
        let mut status: CplexStatus = 0;
        let lp_name: CString = CString::new(name).unwrap();
        let inner = unsafe {
            let ptr = cplex::CPXcreateprob(env.inner, addr_of_mut!(status), lp_name.as_ptr());
            if ptr.is_null() {
                return Err(CplexError::new_ptr(env.inner, status));
            }
            ptr
        };
        Ok(Self {
            // FIXME: Argh, aliasing the ptrs.
            env,
            inner,
        })
    }

    pub fn num_rows(&self) -> CplexInt {
        unsafe { cplex::CPXgetnumrows(self.env.inner.cast_const(), self.inner) }
    }

    pub fn num_cols(&self) -> CplexInt {
        unsafe { cplex::CPXgetnumcols(self.env.inner.cast_const(), self.inner) }
    }

    pub fn set_objective_sense(&mut self, max: bool) -> Result<(), CplexError> {
        let sense = if max {
            cplex::CPX_MAX as CplexInt
        } else {
            cplex::CPX_MIN as CplexInt
        };
        let status = unsafe { cplex::CPXchgobjsen(self.env.inner.cast_const(), self.inner, sense) };
        if status != 0 {
            return Err(CplexError::new(&self.env, status));
        }
        Ok(())
    }

    pub fn new_cols(&mut self, cols: &[CplexCol]) -> Result<(), CplexError> {
        let obj = cols.iter().map(|col| col.obj).collect::<Vec<_>>();
        let lb = cols.iter().map(|col| col.lb).collect::<Vec<_>>();
        let ub = cols.iter().map(|col| col.ub).collect::<Vec<_>>();
        // If all the xctypes are continuous, then you can pass null.
        // "If types is specified, the problem type will be a MIP, even if all variables are specified to be continuous."
        // So double check that if you want a continuous problem, you pass None.
        let xctype = if let None = cols.iter().find(|col| col.xctype != Xctype::Continuous) {
            None
        } else {
            Some(cols.iter().map(|col| col.xctype as i8).collect::<Vec<_>>())
        };

        // Cplex wants a *mut *mut c_char, so we copy the strings into a locally owned buffer.
        let mut col_names = cols
            .iter()
            .map(|col| {
                let mut bytes = col.name.as_bytes().to_vec();
                bytes.push(0);
                bytes
            })
            .collect::<Vec<_>>();
        let mut name_ptrs = col_names
            .iter_mut()
            .map(|name| name.as_mut_ptr().cast::<i8>())
            .collect::<Vec<_>>();
        let status = unsafe {
            cplex::CPXnewcols(
                self.env.inner.cast_const(),
                self.inner,
                cols.len() as CplexInt,
                obj.as_ptr(),
                lb.as_ptr(),
                ub.as_ptr(),
                xctype.map_or(std::ptr::null(), |x| x.as_ptr()),
                name_ptrs.as_mut_ptr(), // Seems weird this is a *mut ptr rather than *const ptr, but whatever.
            )
        };
        if status != 0 {
            return Err(CplexError::new(&self.env, status));
        }
        Ok(())
    }

    /// Add rows to the problem.
    /// Note that it is expected that new_cols has been called before this.
    /// Adding new columns by adding rows is deprecated.
    pub fn add_rows(&mut self, rows: &[CplexRow]) -> Result<(), CplexError> {
        // Don't add new columns here
        let ccnt = 0;
        let rcnt = rows.len() as CplexInt;

        let rhs = rows.iter().map(|row| row.sense.rhs()).collect::<Vec<_>>();
        let sense = rows
            .iter()
            .map(|row| row.sense.as_char())
            .collect::<Vec<_>>();
        // Cplex wants a *mut *mut c_char, so we copy the strings into a locally owned buffer.
        let mut row_names = rows
            .iter()
            .map(|row| {
                let mut bytes = row.name.as_bytes().to_vec();
                bytes.push(0);
                bytes
            })
            .collect::<Vec<_>>();
        let mut name_ptrs = row_names
            .iter_mut()
            .map(|name| name.as_mut_ptr().cast::<i8>())
            .collect::<Vec<_>>();

        let mut rmatbeg = Vec::new();
        let mut rmatind = Vec::new();
        let mut rmatval = Vec::new();

        let mut range_indices = Vec::new();
        let mut range_values = Vec::new();

        for (i, row) in rows.into_iter().enumerate() {
            rmatbeg.push(rmatind.len() as CplexInt);
            for (col, val) in &row.coeffs {
                rmatind.push(*col);
                rmatval.push(*val);
            }
            if row.sense.range().is_some() {
                range_indices.push(i as CplexInt);
                range_values.push(row.sense.range().unwrap());
            }
        }

        let status = unsafe {
            cplex::CPXaddrows(
                self.env.inner.cast_const(),
                self.inner,
                ccnt,
                rcnt,
                rmatval.len() as CplexInt,
                rhs.as_ptr(),
                sense.as_ptr(),
                rmatbeg.as_ptr(),
                rmatind.as_ptr(),
                rmatval.as_ptr(),
                std::ptr::null_mut(),   // We don't add columns here.
                name_ptrs.as_mut_ptr(), // But why does this need to be a *mut ptr?
            )
        };
        if status != 0 {
            return Err(CplexError::new(&self.env, status));
        }

        if !range_indices.is_empty() {
            let status = unsafe {
                cplex::CPXchgrngval(
                    self.env.inner.cast_const(),
                    self.inner,
                    range_indices.len() as CplexInt,
                    range_indices.as_ptr(),
                    range_values.as_ptr(),
                )
            };
            if status != 0 {
                return Err(CplexError::new(&self.env, status));
            }
        }
        Ok(())
    }

    pub fn opt(&mut self) -> Result<(), CplexError> {
        let status = unsafe { cplex::CPXlpopt(self.env.inner.cast_const(), self.inner) };
        if status != 0 {
            return Err(CplexError::new(&self.env, status));
        }
        Ok(())
    }

    pub fn get_solution(&self) -> Result<CplexSolution, CplexError> {
        let n_cols = self.num_cols();
        let n_rows = self.num_rows();

        let mut solstat: CplexInt = 0;
        let mut objval: CplexDouble = 0.0;

        let mut x: Vec<CplexDouble> = vec![0.0; n_cols as usize];
        let mut pi = vec![0.0; n_rows as usize];
        let mut slack = vec![0.0; n_rows as usize];
        let mut dj = vec![0.0; n_cols as usize];

        let status = unsafe {
            cplex::CPXsolution(
                self.env.inner.cast_const(),
                self.inner.cast_const(),
                addr_of_mut!(solstat),
                addr_of_mut!(objval),
                x.as_mut_ptr(),
                pi.as_mut_ptr(),
                slack.as_mut_ptr(),
                dj.as_mut_ptr(),
            )
        };
        if status != 0 {
            return Err(CplexError::new(&self.env, status));
        }
        Ok(CplexSolution {
            solstat,
            objval,
            x,
            pi,
            slack,
            dj,
        })
    }
}

impl<'a> Drop for CplexLp<'a> {
    fn drop(&mut self) {
        println!("Dropping CplexLp");
        if self.inner.is_null() {
            return;
        }
        unsafe {
            let status = cplex::CPXfreeprob(self.env.inner.cast_const(), addr_of_mut!(self.inner));
            if status != 0 {
                let err = CplexError::new(&self.env, status);
                eprintln!("CPXfreeprob failed with error {err:?}");
            }
        }
    }
}

impl CplexError {
    pub fn new(env: &CplexEnv, status: CplexStatus) -> Self {
        Self::new_ptr(env.inner.cast_const(), status)
    }

    pub fn new_ptr(env: cplex::CPXCENVptr, status: CplexStatus) -> Self {
        let mut msg = [0u8; cplex::CPXMESSAGEBUFSIZE as usize];
        unsafe {
            cplex::CPXgeterrorstring(env, status, msg.as_mut_ptr().cast());
        }
        let msg_len = msg.iter().position(|&c| c == 0).unwrap_or(msg.len());
        let msg = String::from_utf8_lossy(&msg[..msg_len]).to_string();
        Self { status, msg }
    }
}
