AWS Tutorial
============

.. contents:: Contents
   :local:

Creating an AWS Account
***********************

The first step is to create an Amazon AWS account if you don’t have one
already.

Go to https://aws.amazon.com and look for an option to create an
account.

AWS Basics
**********

Services
--------

AWS is divided up into many services.

Most of these are enterprise-grade features that we don’t have to worry
about. The two services you should know about are:

-  **EC2**: “Elastic Compute Cloud”, this is the service for renting
   compute time on the servers in Amazon’s cloud
-  **S3**: “Simple Storage Service”, this service is for storing data

Everything we’ll be doing here is within the **EC2** service.

To get to the EC2 dashboard, after logging in:

-  Click ‘Services’ on the toolbar in the upper-left
-  Navigate to “EC2” under the “Compute” section

Regions
-------

Most things you work with in AWS will be scoped to a ‘Region’.
Basically, Amazon has datacenters all over the world and you can choose
which one you want to use for renting servers and storing data.

If you were running a large website, you’d want servers in many regions
for speed and redundancy.

In our case, we just need to be aware that regions exist and how to
select the ‘current’ region.

Whenever you view a service, you are only showed the servers and data
that are in the ‘current’ region. You can see your region and switch it
by looking at the dropdown in the upper-right corner (shown below).

The choice of region isn’t permanent in any sense - you can create
servers in one region, then create some more in another region, and
switch back and forth as you manage either group. However, best to just
pick a region close to you and remain in it. If you ever sign on and all
your servers seem to be missing - don’t fret, you probably just have the
wrong region selected.

Instances and Images
--------------------

To understand why we need to create an Image, it’s important to
understand the AWS terminology:

-  **Instances**: Each server (cloud computer) you provision (rent) is
   called an *instance*.
-  **Volume**: A *volume* is the storage that can be attached to an
   instance. Think of it like plugging a hard drive into a computer.
   Every instance needs at least one volume that contains its operating
   system.
-  **Image**: An *image* is a saved volume that can be used to launch an
   instance.

When you launch an instance, you select an image. This image is copied
onto the volume that is attached to the instance and then the instance
is booted up. The same image can be used to launch multiple volumes.

You don’t have to create your own image - Amazon has many images with
different operating systems already prepared. However, if you use one of
these, you’ll have to install any extra software every time you create
an instance. Instead, we’re going to create an Image with Compass and
its dependencies installed on it so that we only have to do this once.

Note: Images are also referred to as **AMI**\ s (Amazon Machine Images).

Setting up the Compass Image (only need to do this once)
********************************************************

We are going to create an Image with Compass and its dependencies
installed.

To do this, first we are going to create a server with a basic Ubuntu
instance. Install the software. Then save it as an Image.

Creating an Instance
--------------------

-  Navigate to EC2
-  Select ‘Instances’ in the left menu-bar
-  Click the ‘Launch Instance’ button in them middle

Here we choose the base AMI for our instance. Scroll until you see
Ubuntu Server 18.04 and click ‘Select’ on the right.

Next we choose an instance type. Here we get to decide how powerful our
machine is. The caveat is that more powerful machines cost more per
hour. To see pricing, follow `this
link <https://aws.amazon.com/ec2/pricing/on-demand/>`__. Since we are
just installing software, let’s choose a lower-end instance, the
‘t2.large’. Then click ‘Next: Configure Instance Details’ in the bottom
right.

There are a lot of options on this page, but you can ignore most of
them. The one that is good to know about is the ‘Request Spot Instances’
option towards the top. Don’t click this now, but in the future, when
running long jobs, you may want to select this option as spot instances
can save a lot of money. See the `appendix item on Spot
Instances <#spot-instances>`__ for info. For now, just click ‘Next’ on
the bottom right.

On this page, you can set the storage space for your instance. For now,
let’s just leave this one at the default of 8GB since we’re only
installing software. Click ‘Next’ on the bottom right.

On this page you can add tags. This is only useful if you have many
servers and you want to organize them all using tags. Click Next.

On this page you can define which ports are open for your instance. By
default, 22 will be open for SSH. There will be a warning that any IP
address can access your instance. If you’d like you can fix this by
specifying your devices IP address on this page to restrict access to
your machine, but this isn’t necessary. Click ‘Review and Launch’ in the
bottom right.

On this summary page click ‘Launch’. Here it asks you to create an
encrypted key pair for SSH. Once you make a key pair, you can use the
same one for future instances. Since this is our first instance, select
‘Create a new key pair’ from the first dropdown, then give your key pair
a name, e.g. ‘AWSKey’, and click “Download Key Pair”. Make sure to save
the Key Pair .pem file somewhere where you won’t lose it. You can only
download the key pair once, but you can always create new key pairs in
the future if you lose the file.

Now click ‘Launch Instances’. On the next page, you can click ‘View
Instances’ on the bottom to go back to the EC2 >> Instances page where
you can see the status of your instance as it’s booting.

Logging into your Instance
--------------------------

On the EC2 >> Instances page, you can see your instance. To connect to
it via SSH, right-click the instance and click ‘Connect’. You will then
be shown a popup with instructions on how to use your key pair to SSH
into the instance. In this case, I am using this command to connect:

::

   ssh -i "~/AmazonKeys/AWSKey.pem" ubuntu@ec2-35-167-139-94.us-west-2.compute.amazonaws.com

Note that you must provide your own path to the key pair (AWSKey above)
you created. You’ll also have to make sure it’s permissions are ‘400’
(AWS provides instructions on this too on the same popup).

Installing Python and Python Dependencies
-----------------------------------------

The bare Ubuntu 18.04 instance we launched has Python 3.6 installed
already, but we’ll need to install ‘pip’ to download other packages:

::

   sudo apt-get update
   sudo apt-get install python3-pip

Now, let’s install the Compass python package along with all Python
dependencies:

::

   pip3 install numpy
   pip3 install git+https://github.com/yoseflab/Compass.git
   source ~/.profile

That last line just addes the ~/.local/bin directory to the path so
compass can be run more easily.

Installing CPLEX
----------------

CPLEX is an optimization engine that Compass uses and it must be
installed from IBM.

First navigate to the main page for CPLEX
`here <https://www.ibm.com/products/ilog-cplex-optimization-studio>`__

Compass requires the full edition of CPLEX. However, if you are an
Academic (Student or Faculty) you can get a copy for free.

Click the ‘Get student and faculty editions’ link:

Register an account or Log In to an existing account. Once logged in,
scroll down to the ‘Software’ section and to a download page for “ILOG
CPLEX Optimization Studio” v12.9.

Here there are a bunch of options. Find the one named “IBM ILOG CPLEX
Optimization Studio 12.9 for Linux x86-64 Multilingual (CNZM2ML)” and
check it. Then click ‘Download Now’ at the bottom of the page. Then you
probably will need to click the ‘install / re-install Download Director’
popup on the bottom and follow the instructions there.

Finally, you will be able to download
``cplex_studio129.linux-x86-64.bin``. Now we have to send this file to
our running server. To do this, we’ll use scp. One *your computer* (not
the server) run this command, filling in your own paths to files and
your own ubuntu@e2-… command that you used when SSHing into the server.
Notice that we add the ``:~`` designation at the end to send this file
to the home directory on your server instance.

::

   scp -i "<path-to-key>/AWSKey.pem" <path-to-cplex-bin-file> ubuntu@ec2-35-167-139-94.us-west-2.compute.amazonaws.com:~

Once the file has been transferred, on the Ubuntu instance, execute the
installer by running:

::

   sudo apt-get install default-jre

To first install Java (required by CPLEX). Then:

::

   cd ~
   chmod +x cplex_studio129.linux-x86-64.bin
   sudo ./cplex_studio129.linux-x86-64.bin

Follow the instructions in the installer, accepting the license
agreement and choosing to install to the default path
``/opt/ibm/ILOG/CPLEX_Studio129``.

Afterwards, if it has installed succesffully, remove the installer file
with ``rm ~/cplex_studio129.linux-x86-64.bin``

Lastly, we need to install the Python module that comes with cplex. To
do this, run these commands:

::

   cd /opt/ibm/ILOG/CPLEX_Studio129/cplex/python/3.6/x86-64_linux
   sudo python3 setup.py install

Testing everyting out
---------------------

Now to test if everything is installed, simply run:

::

   compass -h

And you should see the help text print out.

Saving the image
----------------

Now that we have everything set up, we’re going to save our instance as
an AMI. This will let us launch instances in the future with all the
software already installed.

On the EC2 Console (in your web browser), right click you instance and
select “Create Image”

Give the image a name and description and hit ‘Create Image’.

Now, on the left side of the dashboard, you can select ‘AMIs’ under the
‘Images’ heading and see the status of your image being created. When
the status goes from ‘pending’ to ‘available’, then you can launch new
instances from this Image.

Running an Analysis
*******************

Now that we have a Compass Image, we can run an analysis in the Amazon
cloud

Launching a (more powerful) instance with our AMI
-------------------------------------------------

Let’s launch an instance using our AMI. First navigate to the EC2 >>
Instances page and click the ‘Launch Instance’ button at the top.

Instead of choosing the Ubuntu 18.04 AMI, we are going to select our
own. Click the MY AMIs tab on the left and select the AMI we just
created by clicking ‘select’ next to it on the right.

On the next page, instead of selecting the low-pwered t2.large instance
type, we are going to select something with much more processing power.
Scroll down and select the c5.24xlarge instance.

Click through the next page until you get the to ‘Add Storage’ page.
Here, up the storage a bit higher so we have more room for data - 20 GB
should be fine.

Raising utilization limit for instance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You *may* get an error launching the instance - something about your
utilization limit for that instance being exceeded. This happens because
for some instance types, Amazon automatically sets a limit of 0 intances
on new AWS accounts. I believe this is done for security purposes (so
that someone getting your login key can’t just spin up 1000 instances to
mine bitcoin at your expense). If you run into this error, it’s fairly
easy to request a limit increase in the support center. I’ve found that
usually Amazon responds to these requests in a few hours.

Click ‘Support’ in the upper-right corner and select ‘Suppert Center’
from the dropdown.

Create a support case and select ‘Service limit increase’. Here select
the Region you are using and the Instance tpe (c5.24xlarge) and request
to raise the instance limit to 1.

Sending data to the instance
----------------------------

Once your instance is running, before logging into it, we are going to
send our expression matrix to the instance using ``scp``.

Remember, you can view connection details for your instance by
right-clicking it and selecting ‘Connect’.

If your ssh command looks like this:

::

   ssh -i "~/AmazonKeys/AWSKey.pem" ubuntu@ec2-54-212-82-196.us-west-2.compute.amazonaws.com

Then your scp command should look like this (assuming your expression
matrix is in a file called data.txt)

::

   scp -i "~/AmazonKeys/AWSKey.pem" data.txt ubuntu@ec2-54-212-82-196.us-west-2.compute.amazonaws.com:~

You may get an error “please log in as ‘ubuntu’ instead of ‘root’”. If
this happens, replace ‘root’ in your scp/ssh command with ‘ubuntu’ and
repeat.

Running the analysis
--------------------

To run the analysis, first log into the instance using SSH.

Next, we are going to start a ``tmux`` session on our instance so that
we can log off of it while the algorithm continues to run.

To start a tmux session, type ``tmux`` at the BASH prompt and hit enter.

Now, follow the instructions on our `Github
page <https://github.com/yoseflab/compass>`__ to start the analysis. Use
the ‘Running on a Single Machine’ option and select 48 processes (half
the number of vCPUs is optimal).

The algorithm will run approximately 48 cells per hour, so if you don’t
want to wait around, you can simply close your terminal window here.
Since we are running inside a ``tmux`` session, the process will
continue.

When logging back in later, to re-attach the ``tmux`` session, simply
run ``tmux attach``.

Downloading the results
-----------------------

When Compass has completed, the outputs are stored in ``reactions.txt``,
``secretions.txt`` and ``uptake.txt``. These should be in the directory
where compass was run by default (can be overridden with the
``--output-dir`` option), which is the home directory in this tutorial.

Now, before we shut down our server, we need to download the results. On
*your computer*, use scp to do this. The command should look like this:

::

   scp -i "~/AmazonKeys/AWSKey.pem" ubuntu@ec2-54-212-82-196.us-west-2.compute.amazonaws.com:~/reactions.txt .

Again, you’ll need to fill in the correct path to your AWSKey.pem file
and your particular instance address. Repeat this two more times for the
``secretions.txt`` and ``uptake.txt`` files.

Cleanup
-------

Now we need to shut down our instance. You will be charged as long as
the instance is in the ‘running’ state, so this is VERY IMPORTANT.

On the EC2 >> Instances page, right click your instance, select
‘Instance State >> Terminate’. You can also ‘Stop’ an instance (unless
it’s a Spot instance) which will let you ‘Start’ it again later without
incurring compute charges while it’s stopped.

Note: You’ll continue to pay a small monthly fee to store the AMI for
Compass we created (about $0.40/month). If you know you won’t need it
again, you can delete it on the EC2 >> AMIs page by right clicking the
AMI and selecting ‘Deregister’

Appendix
********

Spot Instances
--------------

Spot instances are available for much cheaper (often around 1/3 the
cost, `pricing here <https://aws.amazon.com/ec2/spot/pricing/>`__) but
they come with two caveats. First, the price of spot instances can
increase with demand. When you request spot instances you specify what
you are willing to pay and if the price goes above that, then your
instances will be halted. However, if you just bid the current on-demand
price (the default if you leave the bid price blank), then this is not
likely to happen. The second caveat is that you can’t Stop and Start a
spot instance. With on-demand instances you can ‘Stop’ the instance,
essentially pausing it indefinitely. You won’t be charged for compute
costs while it is stopped and then you can start it again at a later
time. With a spot instance, once you stop it you can’t start it again,
so you must download the results of a computation before stopping the
instance. Alternately, you can attach an extra storage volume to the
instance and save your results on that volume, then mount that with
another instance later to get a similar kind of behavior - just with a
bit more work.
