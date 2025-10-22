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

 - Clustering - skip that for now
 - GSMM - let's do that.
    - Recon2 or Recon3.
    - Resolving Gene Names. Not fun,
 - Construct optimization problem.
    - For now, do CPLEX? Perhaps make a trait so that I can drop in whatever backend.

### GEM (GSMM)
 - Why is the abberviation for Genome Scale Metabolic Model = GEM rather than GSMM?
 - We want mouse model, yes? Recon2 appears to be a human one https://www.nature.com/articles/nbt.2488
 - Why does http://humanmetabolism.org/, the link in the Recon2 paper just show a straight up nginx proxy?
 - Formats?
    - _mat ie matlab? Mostly it seems to be our format devised with jsons.
    - _xml ie SBML, a particular schema of XML.
 - Which one to use?
    - I should just dredge up CLI params from somewhere
    - Looks like Recon2, so go with that one.
    - We'll use SBML as that seems to be the one most models are published in.
    - It appears that this crate https://github.com/carrascomj/rust_sbml will suffice.
    - Maybe, perhaps I can use my own XML code to figure it out. 
    - XML crates:
        - quick-xml (used by sbml crate). Per benchmarks on roxmltree, it's probably faster.
            - The unofficial Azure REST API uses it? So I guess maybe just use SBML crate.
        - xmlparser - same author as memmap2. Mostly used via other tools, like aws-smithy-xml or roxmltree
        - roxmltree - read only xml tree. A bit slower than quick-xml, but you get the whole doc.
    - Where is the Recon2 SBML? Hmm, I see Recon2.2 xml gz.
        - Oh the rust_sbml code is actually quite short? Just serde + quick-xml.


### Models
- Recon1
    - rust_sbml parses this one
- Recon2.2
    - rust_sbml chokes on this one due to 'dc:creator'
    - Possibly due to the http://www.sbml.org/sbml/level2/version4 rather than level 2 - version 2
- Recon2_mat
    - What model is this actually?. I kind of suspect it is the same as Recon2.2.
    - Maybe not, I just grepped for R_3HPVSCOAitx and found nothing.
- Perhaps just write code for whatever GSMM I can use for now?


We'll use Recon1 for now perhaps. At least I have a _mat and an xml format that work. Note that rust_sbml Model and ModelRaw have minimal performance differences, so just use Model. The number of reactions and species/metabolites checks out, but this one does not include genes at all, so I may be reduced to parsing things myself. There is a list of associations for a reaction? Also see https://github.com/carrascomj/rust_sbml/issues/2. Yeah it dooes not appear to support gene product associations.


### Parsing MAT
Hmm, the rules part looks like a pain, mostly because the format is not clear to me. In particular, why certain tokens are (x([0-9]+)) or just x([0-9]+) and what's the precedence for & or |. I suppose if it's just left to right precedence, that is okay.

The lack of associativity of the gene rules is a problem I think. Or I guess I should say it's that taking the mean is not associative. I may just try parsing the SBML instead, it has an unambiguous grammar. I suppose an alternative is to collapse the AST so that multi-ORs can be collapsed.

Seems to mostly work though, at least it produces valid stuff I think.

This still seems suspect to me, I'd prefer to see if it's CNF vs DNF. Or perhaps just do what the python code does to check.

Ok, so checking at least the reaction 2OXOADOXm in recon 1, we can check the XML sourced right from SBML against the results of the parsers for the mat format. The compass parsing of the mat format results in what appears to be a mathcing, albeit deeper than neccesary tree, while my shunting-yard esque algorithm results in a different precedence than the SBML. It might be possible to adjust the shunting yard to just prefer ORs to the ANDs. A quick scan of the xml makes it seem like it's all DNF, ie OR of ANDs.

### Python stuff
Using miniforge'd conda. install:
```
 numpy pandas python-libsbml
```

## Revisiting this
So initial thoughts are that I should just stop using the binary tree stuff? Just use the arbitrarily shaped tree. Seems like that will work better. After double checking the SBML, I can see very clearly things like this in Recon1:
``` xml
        <fbc:geneProductAssociation xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2">
          <fbc:or sboTerm="SBO:0000174">
            <fbc:geneProductRef fbc:geneProduct="G_2572_AT1" />
            <fbc:geneProductRef fbc:geneProduct="G_51380_AT1" />
            <fbc:geneProductRef fbc:geneProduct="G_2571_AT1" />
            <fbc:geneProductRef fbc:geneProduct="G_2571_AT2" />
          </fbc:or>
        </fbc:geneProductAssociation>
```
## Setting up my venv again
### Compass
From the repo root
```bash
conda create -n compass_env python=3.12 
conda activate compass_env
conda install numpy pandas python-libsbml
pip install .
```
### Cplex
Hmm, it seems like the community version lacks the python directory? That will be a problem. Oh, you can install that one with pip it looks like. So with the community version (which probably won't work with compass anyways)
```
conda install ibmdecisionoptimization::cplex
```
Then you can at least do `python -m compass.main`

### cuOPT
Note you need cuda installed. 

I am simply following the [docs](https://docs.nvidia.com/cuopt/user-guide/latest/cuopt-c/quick-start.html)
```
# CUDA 13 - C library
conda install -c rapidsai -c conda-forge -c nvidia libcuopt=25.10.* cuda-version=13.0
# CUDA 13 - python
conda install -c rapidsai -c conda-forge -c nvidia cuopt=25.10.* cuda-version=13.0
``` 
Then I find the C headers with a
```sh
find $CONDA_PREFIX/include -type f -name '*.h' | grep cuopt
```
For replicability, I have this package
```sh
$ conda list cuopt
# packages in environment at /home/user/miniforge3/envs/compass_env:
#
# Name                    Version                   Build  Channel
cuopt                     25.10.00        cuda13_py312_251014_99e549ce    nvidia
cuopt-mps-parser          25.10.00        py312_251014_99e549ce    nvidia
libcuopt                  25.10.00        cuda13_251014_99e549ce    nvidia
```

### Getting libcusparse
Then as I am using WSL 2 I will also need to [install CUDA toolkit](https://docs.nvidia.com/cuda/wsl-user-guide/index.html). In this case, go to https://developer.nvidia.com/cuda-downloads and select the appropriate version (in my case: Operating System: Linux; Architecture: x86_64; Distribution: WSL-Ubuntu; Version: 2.0).

For ease of the reader:
```
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-13-0
```

Then the [Post-installation actions](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#post-installation-actions)
```sh
export PATH=${PATH}:/usr/local/cuda-13.0/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/cuda-13.0/lib64
```

## Linking debugging
Usually I try to statically link things to avoid issues with dynamic linking, but given nvidia only distributes dynamic libraries, static linking is not very practical.

1. libcuopt.so: "error while loading shared libraries: libcuopt.so: cannot open shared object file: No such file or directory"
   1. This library is installed at $CONDA_PREFIX/lib. On linux use the `LD_LIBRARY_PATH` (on Mac it should be `DYLD_LIBRARY_PATH`) variable: `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib cargo run`. Note that [conda discourages settings this persistently](https://docs.conda.io/projects/conda-build/en/stable/resources/use-shared-libraries.html#shared-libraries-in-macos-and-linux), as it can interfere with resolving system libraries (ie /usr/lib/), so I generally set the variable only for the commands that require it.
1. libcupsparse.so : error while loading shared libraries: libcusparse.so.12: cannot open shared object file: No such file or director
    1. These are installed with the CUDA toolkit.


# Compass algorithm

1. Get gene expression data - I should have some sitting around
1. Parse - probably into polars data frame (hmm, where did I put that code. I should have pushed that stuff to git)
1. Group gene symbols
1. Isoform summing and dealing with reversible reactions
1. Simple single cell penalty computation
1. Compute max throughput of reactions
1. Construct flux balance analysis problem and solve

Optional steps
1. Cache maximum flux through model
    1. Probably required for even half decent performance
1. Smoothing, over lambda > 0
    1. Latent space input: simple
    1. PCA: I can implement that fine.
    1. tsne: Tricky to do, ideally find a library.
    1. knn: Hard to do it fast, ideally find a library.

## Some design thoughts:

1. Probably want to use python calling into Rust. Which parts should be py and which rs?
    1. CLI? I like clap, but probably should be python, as it's the entry point. Argparse stil suffices.
    1. Data loading - eh we'll use polars in python may as well.
    1. Data preprocessing also python to make use of the various packages
    1. Do gsmm part in rust as a polars extension - should be interesting
    1. Caching doesn't make a big difference
    1. Optimization engine
        1. So obviously writing one myself would be Rust
        1. GLPK has a C interface, while technically C could do it, it's probably easier for me to use Rust so I can ensure the type layouts match.
        1. Hmm, probably call into the gsmm rust module to get numpy arrays to feed to the optimization engines.

### Misc: cuOPT error
I have gotten this error probably 1 in 3 runs.
```sh
terminate called after throwing an instance of 'raft::cusparse_error'
  what():  cuSparse error encountered at: file=/tmp/conda-bld-output/bld/rattler-build_libmps-parser/work/cpp/src/dual_simplex/sparse_matrix_kernels.cuh line=134: call='cusparseSpGEMM_compute(handle->get_cusparse_handle(), CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, cusparse_data.alpha.data(), cusparse_data.matA_descr, cusparse_data.matDAT_descr, cusparse_data.beta.data(), cusparse_data.matADAT_descr, CUDA_R_64F, CUSPARSE_SPGEMM_ALG3, cusparse_data.spgemm_descr, &cusparse_data.buffer_size_2_size, cusparse_data.buffer_size_2.data())', Reason=7:internal error
Obtained 11 stack frames
#1 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: raft::cusparse_error::cusparse_error(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&) +0x5a [0x7ac15f62ce0a]
#2 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so(+0x89c11f) [0x7ac15fc9c11f]
#3 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: cuopt::linear_programming::dual_simplex::iteration_data_t<int, double>::form_adat(bool) +0x4ad [0x7ac15fcc08bd]
#4 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: cuopt::linear_programming::dual_simplex::iteration_data_t<int, double>::iteration_data_t(cuopt::linear_programming::dual_simplex::lp_problem_t<int, double> const&, int, cuopt::linear_programming::dual_simplex::simplex_solver_settings_t<int, double> const&) +0x295d [0x7ac15fcd2d5d]
#5 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: cuopt::linear_programming::dual_simplex::barrier_solver_t<int, double>::solve(double, cuopt::linear_programming::dual_simplex::barrier_solver_settings_t<int, double> const&, cuopt::linear_programming::dual_simplex::lp_solution_t<int, double>&) +0x26c [0x7ac15fcdc22c]
#6 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: cuopt::linear_programming::dual_simplex::lp_status_t cuopt::linear_programming::dual_simplex::solve_linear_program_with_barrier<int, double>(cuopt::linear_programming::dual_simplex::user_problem_t<int, double> const&, cuopt::linear_programming::dual_simplex::simplex_solver_settings_t<int, double> const&, cuopt::linear_programming::dual_simplex::lp_solution_t<int, double>&) +0xad4 [0x7ac15fd3b294]
#7 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: std::tuple<cuopt::linear_programming::dual_simplex::lp_solution_t<int, double>, cuopt::linear_programming::dual_simplex::lp_status_t, double, double, double> cuopt::linear_programming::run_barrier<int, double>(cuopt::linear_programming::dual_simplex::user_problem_t<int, double>&, cuopt::linear_programming::pdlp_solver_settings_t<int, double> const&, cuopt::timer_t const&) +0x4de [0x7ac15f631a9e]
#8 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/cuopt/linear_programming/solver/../../../../../libcuopt.so: void cuopt::linear_programming::run_barrier_thread<int, double>(cuopt::linear_programming::dual_simplex::user_problem_t<int, double>&, cuopt::linear_programming::pdlp_solver_settings_t<int, double> const&, std::unique_ptr<std::tuple<cuopt::linear_programming::dual_simplex::lp_solution_t<int, double>, cuopt::linear_programming::dual_simplex::lp_status_t, double, double, double>, std::default_delete<std::tuple<cuopt::linear_programming::dual_simplex::lp_solution_t<int, double>, cuopt::linear_programming::dual_simplex::lp_status_t, double, double, double> > >&, cuopt::timer_t const&) +0x27 [0x7ac15f631fc7]
#9 in /home/bschel/miniforge3/envs/compass_env/lib/python3.12/site-packages/numpy/_core/../../../../libstdc++.so.6(+0xd828c) [0x7ac30509b28c]
#10 in /lib/x86_64-linux-gnu/libc.so.6(+0x9caa4) [0x7ac30889caa4]
#11 in /lib/x86_64-linux-gnu/libc.so.6(+0x129c6c) [0x7ac308929c6c]
```

Another error:
```
CUDA Error detected. CUDA Error detected. cudaErrorIllegalAddress cudaErrorIllegalAddress an illegal memory access was encounteredan illegal memory access was encountered

python: /tmp/conda-bld-output/bld/rattler-build_libmps-parser/host_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_pl/include/rmm/mr/device/cuda_memory_resource.hpp:80: virtual void rmm::mr::cuda_memory_resource::do_deallocate(void*, std::size_t, rmm::cuda_stream_view): Assertion `status__ == cudaSuccess' failed.
```
But after retrying with CUDA_LAUNCH_BLOCKING=1, it seems to work.

When I use CUDA_LAUNCH_BLOCKING=0, I also get PrimalInfeasible issues as well as crashes. This seems like either I am doing something wrong/unexpected or they are.

The cuda computer-sanitizer tool seems to just be too slow, no problem completed after about an hour, though it does at least start the solver. At least with memcheck. synccheck, racecheck, and initcheck all also take a significant amount of time.

### Memory usage
Seems to be some memory usage issues from running it long enough to cache everything in one process. On my laptop, it ended up using a lot more RAM than availible and swapping to disk significantly. Probably the simplest fix is to spawn multiple processes that process a smaller portion. Unsure if this is cuOpt leaking memory or if it's python.