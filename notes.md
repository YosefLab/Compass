# Notes on install

- Need sudo apt install JRE?
- Need sudo for the defautl install path? Perhaps this might fuck with path and env?
- Which python version to use for it?
    - I'll use python3 /opt/ibm/ILOG/CPLEX_Studio221/python/setup.py install
- Missing distutils with system python3 on Ubuntu?
    - sudo apt install python3-distutils
    - /opt/ibm/ILOG/CPLEX_Studio221/cplex/python/3.10/x86-64_linux/setup.py:15: DeprecationWarning: The distutils package is deprecated and slated for removal in Python 3.12. Use setuptools or check PEP 632 for potential alternatives 
    - from distutils.core import setup /usr/lib/python3.10/distutils/dist.py:274: UserWarning: Unknown distribution option: 'zip_safe'
- Also, I needed to use sudo python3 /opt/ibm/ILOG/CPLEX_Studio221/python/setup.py install
    - Else, I got 
```
running install
running build
running build_py
creating build
error: could not create 'build': Permission denied
```
- Sure enough, python3; import cplex; works.

## Compass install
Apparently I also need to `sudo apt install python3-setuptools`
Use `sudo python3 setup.py install` from the Compass repo.
Jesus, it's got a lot of invalid version warnings.

Now `compass -h` works.

## Testing
Lets go grab the Th17 cells again lol. I can probably find it in my email along with results.
 - Hmm, seems relatively difficult to go find it via grepping for *.csv or *.tsv
 - It may be locatable on my giant-ass hard drive. It seems to have a very funny layout, perhaps due to robocopy shenanigians.
 - Oh, even this machine had it. 

## Development plan

 1. Get cplex working.
 2. Get parts of compass working.
 3. Do GSMM integration.
 4. Make my own LP solver.
    4a. Use ndalgebra?
    4b. Use my FPGA. I imagine this can get pretty fast as you can do something like encode the gsmm directly into the FPGA. Perhaps this will exhaust the gates though?

 ### Getting cplex working
 So, see the cplex headers via `sh find /opt/ibm/ILOG/CPLEX_Studio221/cplex -name *.h`. They are all in `/opt/ibm/ILOG/CPLEX_Studio221/cplex/include/ilcplex/`

 We can find the lib in `/opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/`, we shall ignore the jar because I have no interest in java. See `/opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/libcplex.a`, probably.

 As per usual with bindgen, see https://rust-lang.github.io/rust-bindgen/requirements.html and install libclang-devel.

 Oh yeah I should go check if the cplex-sys crate does what I am doing already. Not really, it seems to run afoul of the same issue I have where FP constants are double-defined. Possibly due to enum colliding with constant? See https://github.com/rust-lang/rust-bindgen/issues/687#issuecomment-316983630
 ``` rs
pub const FP_ZERO: u32 = 2;
pub const FP_ZERO: _bindgen_ty_1 = 2;
 ```
 So I do indeed get to remake the wheel in my very own image.

 So I have the basic LP example working and ranges working, that will do for now.

 ### Now to get parts of compass working.