# Installation and setup

## 1. Requirements
You need to have Python 3.10 or newer installed on your system. If you don't have
Python installed, we recommend installing [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

## 2. Install Compass

```bash
pip install git+https://github.com/YosefLab/Compass.git@main
```

## 3. Obtain a Gurobi WLS license

At the heart of the Compass algorithm is linear programming. Compass computes a penalty score for each reaction in each cell by solving a constrained optimization problem using the Gurobi linear solver. In order to use Compass, you need to obtain a Gurobi WLS license, which is free for academic use.

Please refer to [this link][link-gurobi] to obtain a Gurobi WLS license. If you follow the instructions correctly, you should be able to obtain a `gurobi.lic` file that contains the access ID, secret, and license ID of your WLS license. Please store this file as it is required to set up Compass.

## 4. Setting up your Gurobi WLS license

In order to use Compass, you must provide your Gurobi WLS license to the Gurobi API. To do this, you can run the command
```bash
compass --set-license <PATH_TO_LICENSE>
```
where `<PATH_TO_LICENSE>` is the path to the `gurobi.lic` file you obtained in the previous step. This stores your license information in a default location within Compass. After you have successfully done so, you can proceed to run Compass without the `--set-license` parameter as Compass will directly read your license information from the default location.


[link-gurobi]: https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-license
[link-api]: https://compass_v2.readthedocs.io/latest/api.html