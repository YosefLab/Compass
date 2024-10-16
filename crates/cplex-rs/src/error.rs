use cplex_bindings as cplex;

use crate::{env::CplexEnv, types::CplexStatus};

#[derive(Debug)]
pub struct CplexError {
    #[expect(dead_code, reason = "Used for debugging")]
    status: CplexStatus,
    #[expect(dead_code, reason = "Used for debugging")]
    msg: String,
}

impl CplexError {
    pub fn new(env: &CplexEnv, status: CplexStatus) -> Self {
        Self::new_ptr(env.inner(), status)
    }

    pub(crate) fn new_ptr(env: cplex::CPXCENVptr, status: CplexStatus) -> Self {
        let mut msg = [0u8; cplex::CPXMESSAGEBUFSIZE as usize];
        unsafe {
            cplex::CPXgeterrorstring(env, status, msg.as_mut_ptr().cast());
        }
        let msg_len = msg.iter().position(|&c| c == 0).unwrap_or(msg.len());
        let msg = String::from_utf8_lossy(&msg[..msg_len]).to_string();
        Self { status, msg }
    }
}
