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

## Now for the algorithm

Steps: 
0. Generating the cache is a kind of 0th step.
1. microcluster cells? May be tricky to replicate exactly.
2. Compute penalties
3. Run algorithm core

### Microclustering
I'm generally inclined to skip this for now. It's really an additional step, not neccesarily core to what I want to do.

Oh, it's the VISION algorithm. The only really tricky thing is the leiden algorithm, which I am not at all familiar with.

KMeans is inconvenient, but not neccesarily complicated to implement. It may be complicated to do quickly.

### Penalties
Could use knn smoothing, but once again, we can skip that. Knn is not that hard at least. Tsne is pretty complicated, I don't really want to re-implement that for this project. That would be its own project at least.

The core of the penalties is neccesary though.

1. Read the data
2. Preprocess
 a. Aggregate any duplicate gene symbols
3. Setup model see init_model
4. Actually do the evaluation.

### Algorithm - flux balance analysis is the core of the work.
We want to do this.

### Linear algebra library
Faer - 228k downloads. Pure rust implementation. Probably the coolest one. Does have some conversion traits for the following two.
Nalgebra - 19m downloads. Uses lapack with some configuration it appears.
Ndarray - 18m downloads. Can use BLAS with some configuration.

sprs maybe? For sparse stuff.

Another factor - how compatible with polars? Well, polars can convert things to ndarray. Seems like a reasonable start?
Also, burn has some built-in compatability with ndarray.

Well, what operations do I need:
 1. np.log2
 2. Add 1 to all elements of the array
Pretty easy.

Down the line may want: tsne? Umap?

Whatever, back to more serious concerns:
 1. conda install openblas
 2. sudo apt intsall pkg-config libssl-dev - needed for openssl-sys to find the openssl system version, used by openblas-src for something or another I guess.

### Misc nonsense
How to install compass on my conda pip? Seems like a pita, mostly because sudo obliterates the environment and path.
Hmm, very irritating permission error. Here. 

### Penalty computation
So if I use the default for the and/or functions: mean/sum. Can I use a matrix multiplication for this? I think so.

So mean(x1 + x2 + x3) = x1 / 3.0 + x2 / 3.0 + x3 / 3.0.

So you can do something like:

let mut row = vec![0.0; cols];
let mut stack = vec![(root, 1.0)];
while let Some((t, m)) = stack.pop() {
    match t {
        gene => row[gene.id] += m;
        or => {
            stack.push((child, m));
        }
        and => {
            stack.push((child, m / num_children))
        }
    }
}

Note that for other, non-linear functions you can't do this. Also, it's not clear to me how to efficiently do the multiplication. Apache arrow, I suspect, is not quite designed with this approach in mind.

Also, other ops will require another strategy, so lets do the column function thing instead.

### Penalties debugging
It appears most of my results match the python code. 13DAMPPOX_pos vs my 13DAMPPOX differs though. Also I am not splitting pos and neg. And to be fair, pos vs neg makes no difference for penalties, does it? The gene-protein rule should be the same.
 1. Check that py does generate a bunch of identical penalties for pos vs neg
 1. Check why the 13DAMPPOX differs from 13DAMPPOX_pos

Could it be the isoform summing? My rule for 13DAMPPOX is the same it appears, but I note that AOC2 appears twice. Yeah that fixed it.

Okay the difference appears to be the gene symbols? The first 3 match and then everything else screws up. CYP4F12 vs CYP4F14. Oh is it human vs mouse? Yeah that fixed that the debugging it appears. The only difference I see is in the order of symbols? For python, using a set may result in an arbitrary order?

Python {'CYP3A44', 'CYP3A41A', 'CYP3A41B', 'CYP3A11', 'CYP3A16'} vs ["CYP3A41B", "CYP3A16", "CYP3A11", "CYP3A41A", "CYP3A44"]

Hmm, also note that in Gene eval_expression there is the alt_symbol matching option, which I don't support currently. Mostly because polars is pretty restrictive in comparison and requires nonzero effort vs pandas.

Omg is it just the capitalization? CYP1A1 vs Cyp1a1 lmao. But no, there is already expression.index = expression.index.astype('str').str.upper()  # Gene names to upper. So why is Cyp1a1 not being accepted by the python code?

Hmm, for the same reaction the python code is gettign a different number of genes to scan? 26 vs 15? Oh, my code is not trying to resolve ones without any symbols. I get different values than pandas though. Okay yeah my code is wrong somewhere. I get CYP2C29 with 942.02 rather than 812.79 (python and manually examining the file agree). Oh it looks like: (CYP2C29 + CYP2C38 + CYP2C39) 812.79 +  116.94 + 12.29 = 942.02, so it's because I am including the alt symbols.

So once again consulting the python code, we see that Gene eval_expression will match by name and only use alt symbols if an exact match for the name could not be found, and takes the average across them. While my code just sums across all. Hmm, using polars I can probably change the filtering a bit. Not sure how efficient it will be though.

So I think these new changes should work.

Okay, so it appears my secondary code is not quite cooperating. We have python:
```
ENO1B found alt symbol ENO1 expression 738.88
ENO3 found in index expression 3.43
ENO2 found in index expression 0.45
```
Rust:
```
2025-02-13T13:46:01.116242Z DEBUG compass::penalties: compass/src/penalties.rs:21: Gene GeneId { id: 743 } has expr 0
2025-02-13T13:46:01.116274Z DEBUG compass::penalties: compass/src/penalties.rs:21: Gene GeneId { id: 744 } has expr 3.4299999999999997
2025-02-13T13:46:01.116293Z DEBUG compass::penalties: compass/src/penalties.rs:21: Gene GeneId { id: 746 } has expr 0.45
```
So this is not resolving the alt symbol ENO1

## Plotting?

Some options for plotting things
 1. ggplot - but then I have to use R. I don't want to do that, not for at least 10 years.
 1. seaborn - matplotlib, but in a very, very nice hat. Can do dataframe interchange protocol with polars? I'm also pretty used to seaborn. It's fairly nice.
 1. plotnine - python-ish version of ggplot?
 1. altair - another grammar of graphics thing? It works with polars very nicely. That seems good enough for me.
 1. Plotly - very interactive thing?