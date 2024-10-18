use std::{ffi::CString, ptr::addr_of_mut};

use cplex_bindings as cplex;

use crate::{
    env::CplexEnv,
    error::CplexError,
    types::{CplexCol, CplexDouble, CplexInt, CplexRow, CplexSolution, CplexStatus, Xctype},
};

pub struct CplexLp<'a> {
    inner: cplex::CPXLPptr,
    // TODO: Arc it?
    env: &'a CplexEnv,
}

impl<'a> Drop for CplexLp<'a> {
    fn drop(&mut self) {
        println!("Dropping CplexLp");
        if self.inner.is_null() {
            return;
        }
        unsafe {
            let status = cplex::CPXfreeprob(self.env.inner(), addr_of_mut!(self.inner));
            if status != 0 {
                let err = CplexError::new(&self.env, status);
                eprintln!("CPXfreeprob failed with error {err:?}");
            }
        }
    }
}

impl<'a> CplexLp<'a> {
    pub fn new(env: &'a CplexEnv, name: &str) -> Result<Self, CplexError> {
        let mut status: CplexStatus = 0;
        let lp_name: CString = CString::new(name).unwrap();
        let inner = unsafe {
            let ptr = cplex::CPXcreateprob(env.inner(), addr_of_mut!(status), lp_name.as_ptr());
            if ptr.is_null() {
                return Err(CplexError::new_ptr(env.inner(), status));
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
        unsafe { cplex::CPXgetnumrows(self.env.inner(), self.inner) }
    }

    pub fn num_cols(&self) -> CplexInt {
        unsafe { cplex::CPXgetnumcols(self.env.inner(), self.inner) }
    }

    pub fn set_objective_sense(&mut self, max: bool) -> Result<(), CplexError> {
        let sense = if max {
            cplex::CPX_MAX as CplexInt
        } else {
            cplex::CPX_MIN as CplexInt
        };
        let status = unsafe { cplex::CPXchgobjsen(self.env.inner(), self.inner, sense) };
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
                self.env.inner(),
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
                self.env.inner(),
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
                    self.env.inner(),
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
        let status = unsafe { cplex::CPXlpopt(self.env.inner(), self.inner) };
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
                self.env.inner(),
                self.inner,
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
