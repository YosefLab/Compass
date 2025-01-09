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

### Python stuff
Using miniforge'd conda. install:
```
 numpy pandas python-libsbml
```