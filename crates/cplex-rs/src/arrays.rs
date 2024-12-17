//! How do I want to handle the array of struct -> struct of arrays conversion?
//! A couple of ways to handle this:
//!

struct AStruct {
    a: i32,
    b: i32,
    c: i32,
}

struct SimpleSoa {
    a: Vec<i32>,
    b: Vec<i32>,
    c: Vec<i32>,
}

struct UnsafeSoa {
    /// A pointer to data of format [a_0, ..., a_len, b_0, ..., b_len, c_0, ..., c_len]
    /// Note that there may need to be some padding if a,b,c require different alignments.
    data: *mut i32,
    /// The number of elements
    len: usize,
}

struct UnsafeStridedSoa {
    /// [a_0, a_1, ..., a_n, b_0, b_1, ..., b_n, c_0, c_1, ..., c_n, a_n+1, a_n+2, ..., a_2n, ..., a_len]
    /// Where n is some constant s.t. we have a cache line of a,b,c interleaved.
    data: *mut i32,
    /// The number of elements
    len: usize,
}
