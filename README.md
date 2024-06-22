# Compass

Compass is an algorithm to characterize the metabolic state of cells based on single-cell RNA-Seq and flux balance analysis (FBA). Its development was motivated by the challenge in characterizing the metabolic states of single cells on scale with current metabolic assays. Compass is an in silico approach to infer metabolic status of cells based on transcriptome data, which is readily available even in single cell resolutions.

For instructions on how to install and use Compass, visit the [documentation](https://yoseflab.github.io/Compass/). For a detailed description of the algorithm, refer to Wagner et al., <i>Cell</i> 2021 ([link](https://doi.org/10.1016/j.cell.2021.05.045)).

## Installation and setup

You need to have Python 3.10 or newer installed on your system. If you don't have
Python installed, we recommend installing [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

There are several alternative options to install compass_v2:

1. Install the latest development version:

```bash
pip install git+https://github.com/YosefLab/Compass.git@main
```

## User guide

Please refer to the [documentation][link-docs]. In particular, the

-   [Compass settings documentation][link-settings].

## Release notes

See the [changelog][changelog].

## Contact

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/DPLemonade/compass_v2/issues
[changelog]: https://compass_v2.readthedocs.io/latest/changelog.html
[link-docs]: https://compass-sc.readthedocs.io
[link-repo]: https://github.com/YosefLab/Compass
[link-settings]: https://compass-sc.readthedocs.io/en/latest/settings.html
[link-manuscript]: https://doi.org/10.1016/j.cell.2021.05.045
