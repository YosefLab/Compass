Installing Compass
==================

.. contents:: Contents
   :local:

Requirements
************
 - Python (2.7+ recommended)
 - IBM CPLEX Optimization Studio (12.9+ recommended)

If you are using Python 2.x you will need to install Cplex version 12.9 or earlier as `Cplex 12.10 has dropped support for Python 2.x. <https://www.ibm.com/support/knowledgecenter/SSSA5P_12.10.0/ilog.odms.studio.help/CPLEX/ReleaseNotes/topics/releasenotes12100/convert.html>`__ As of Cplex version 20.1.0 Python 3.7 and 3.8 are supported. See the following sections for more details on installing Cplex.

Older Python or Cplex versions may work but have not been tested.

Installing Compass with Pip
***************************

The simplest way to install Compass is using pip to install from the github repository. This can be done with the following command in the terminal:

.. code:: bash

   pip install git+https://github.com/yoseflab/Compass.git

Now to test if everything is installed, simply run:

.. code:: bash

   compass -h

You should see the help text print out if installation was succesful :) For more details on how to use Compass you can visit our :doc:`tutorial <tutorial>`.


Installing Cplex
****************

CPLEX is an optimization engine that COMPASS uses and it must be
installed from IBM. We use the Python API for Cplex (not Docplex which is a higher level API). 

We note that it may be simpler to download the
installation files locally and SCP/File Transfer them over to the
corresponding server.

Download
--------

First navigate to the main page for CPLEX `here <https://www.ibm.com/products/ilog-cplex-optimization-studio>`__ 

COMPASS requires the full edition of CPLEX, which is free for academic use. Click the ‘Get student and faculty editions’ link.

Register an account or Log In to an existing account. Once logged in, \
scroll down to the ‘Software’ section and to a download page for “ILOG CPLEX Optimization Studio”.

Here there are a bunch of options. Find the one corresponding to your 
operating system, for example “IBM ILOG CPLEX Optimization Studio 12.9
for Linux x86-64 Multilingual (CNZM2ML)” and check it. Then click
‘Download Now’ at the bottom of the page. Then you probably will need to
click the ‘install / re-install Download Director’ popup on the bottom
and follow the instructions there.

Finally, you will be able to download
``cplex_studio129.linux-x86-64.bin`` (or the corresponding file for your version and OS) .

Linux Install Walkthrough
-------------------------

Once the file has been transferred, on the Ubuntu instance, execute the
installer by running the following commands:

First install Java (required by CPLEX) if you haven’t already done so.

.. code:: bash

   sudo apt-get install default-jre

Then:

.. code:: bash

   cd ~
   chmod +x cplex_studio129.linux-x86-64.bin
   sudo ./cplex_studio129.linux-x86-64.bin

Follow the instructions in the installer, accepting the license
agreement and choosing to install to the default path
``/opt/ibm/ILOG/CPLEX_Studio129``.

Afterwards, if it has installed successfully, you can remove the installer file
with 
.. code:: bash

   rm ~/cplex_studio129.linux-x86-64.bin

Lastly, we need to install the Python module that comes with cplex. To
do this, run these commands:

.. code:: bash

   cd /opt/ibm/ILOG/CPLEX_Studio129/cplex/python/3.6/x86-64_linux
   sudo python3 setup.py install

If all is good, cplex will be installed! To test this simply open a
python instance and run the following command

.. code:: bash

   import cplex

If there are no errors, you’re good to go!

Installing on other Operating Systems
-------------------------------------

For installation on Windows or Mac the process will be similar. Navigate to [Cplex Install Directory]/cplex/python/3.6/[OS] and run

.. code:: bash

   python3 setup.py install

For more detailed instructions see `IBM's Knowledge Center <https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html>`__ and 
using the ''Change version or product'' to navigate to the version of Cplex you downloaded. 
Then see `here <https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html>`__ for how to setup the Python API of Cplex.
As before, if you can open a Python instance and run

.. code:: bash

   import cplex

Then you are good to go!