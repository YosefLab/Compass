Installation
==================

.. contents:: Contents
   :local:

Requirements
************
 - python >= 2.7
 - gurobipy >= 11.0.0
 - numpy >= 1.12
 - pandas >= 0.20
 - scikit-learn >= 0.19
 - scipy >= 1.0

Install Compass
*******************
If Numpy is not installed, you can install it with

.. code:: bash

   python -m pip install numpy
   
This needs to be installed before the other requirements because a C extension needs the location of numpy headers to compile.

.. note::

   Accessing pip through python -m pip is done to emphasize that on systems with multiple python installations 
   (e.g. python 2.7 and python 3.6) Compass and Cplex must be installed to the same version of python. 
   It is otherwise identical to just using pip. Using sudo can also invoke a different version of python 
   depending on your environment.

Then simplest way to install Compass is using pip to install from the github repository. This can be done with the following command in the terminal:

.. code:: bash

   python -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

This command will also update Compass to the newest version. Alternatively, you can clone the Compass repository and run setup.py.

Now to test if everything is installed, simply run:

.. code:: bash

   compass -h

You should see the help text print out if installation was succesful. For more details on how to use Compass you can visit our tutorial.

Obtain a Gurobi WLS license
***************************

At the heart of the Compass algorithm is linear programming. Compass computes a penalty score for each reaction in 
each cell by solving a constrained optimization problem using the Gurobi linear solver. 
In order to use Compass, you need to obtain a Gurobi WLS license, which is free for academic use.

Please refer to `this link <https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-license>`__ 
to obtain a Gurobi WLS license. If you follow the instructions correctly, you should be able to obtain a 
``gurobi.lic`` file that contains the access ID, secret, and license ID of your WLS license. 
Please store this file as it is required to set up Compass.

Set up your Gurobi WLS license
******************************

In order to use Compass, you must provide your Gurobi WLS license to the Gurobi API. To do this, you can run the command

.. code:: bash

   compass --set-license <PATH_TO_LICENSE>

where ``<PATH_TO_LICENSE>`` is the path to the ``gurobi.lic`` file you obtained in the previous step. 
This stores your license information in a default location within Compass. 
After you have successfully done so, you can proceed to run Compass without the ``--set-license`` parameter 
as Compass will directly read your license information from the default location.
Refer to the `next section <https://compass-sc.readthedocs.io/en/latest/quickstart.html>`__ of the documentation to learn how to run Compass.