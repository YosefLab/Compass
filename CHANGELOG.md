# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## [Unreleased]

### Added

- Added Module-Compass
- Added [Turbo-Compass](https://compass-sc.readthedocs.io/en/latest/turbo_compass.html)
- Added pseuobulking tutorial
- Added support for Human1 and Mouse1 GEMs
- Added list of core RECON2 reactions for faster computation

### Changed
- Documentation now supported by `Read the Docs`
- Updated cache computation for deterministic results
- Updated postprocessing notebook to use reproducible COMPASS output
- Replaced CPLEX with Gurobi Linear Solver

### Deprecated
- Deprecated microclustering
- Deprecated torque queue