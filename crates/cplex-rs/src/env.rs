use std::ptr::addr_of_mut;

use cplex_bindings::{self as cplex, params::CplexParam};

use crate::{
    error::CplexError,
    types::{CplexInt, CplexStatus},
};

pub struct CplexEnv {
    inner: cplex::CPXENVptr,
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

    pub(crate) fn inner_mut(&mut self) -> cplex::CPXENVptr {
        self.inner
    }

    pub(crate) fn inner(&self) -> cplex::CPXCENVptr {
        self.inner.cast_const()
    }

    pub fn set_param(&mut self, param: CplexParam, value: cplex::CPXINT) -> Result<(), CplexError> {
        let status = unsafe { cplex::CPXsetintparam(self.inner_mut(), param as CplexInt, value) };
        if status != 0 {
            return Err(CplexError::new_ptr(self.inner(), status));
        }
        Ok(())
    }
}

impl Drop for CplexEnv {
    fn drop(&mut self) {
        unsafe {
            let status = cplex::CPXcloseCPLEX(addr_of_mut!(self.inner));
            if status != 0 {
                let err = CplexError::new_ptr(self.inner, status);
                eprintln!("CPXcloseCPLEX failed with error {err:?}");
            }
        }
    }
}
