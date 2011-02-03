This utility can be used to create a rendering cluster on Amazon's EC2
service. Note that so far, it has only only been used on Ubuntu Linux, hence
it might not work on other distributions / operating systems. 

If you have an EC2 account, please look at the beginning of config.py
and fill out the required authentication data. Furthermore, ensure that the
'boto' Python library (named 'python-boto' in Ubuntu) is installed on your
machine, and that all necessary files are located in the same directory
as cluster.py (the SSH key and a .s3cfg file for s3cmd if you are planning 
to use Amazon's S3 service to stream heavy data files)

You can execute

$ ./cluster.py

to get a rough overview of what the program does. The usual approach is to
create a few rendering nodes, e.g.

$ ./cluster.py createNodes c1.xlarge 8

which allocates 8 extra-large high CPU instances. Afterwards, you can install
Mitsuba and all dependencies by running

$ ./cluster.py install

This might generate a few warnings from the installation process, which are
usually safe to ignore. Finally, you can launch Mitsuba on all machines

$ ./cluster.py start

The last message will print the name of a "Head node". This is the machine
which you can register in the Mitsuba GUI (select Tools->Settings->
Network->Plus Button and create a direct connection with port 7554) or on 
the command line, e.g. by running

$ mitsuba -c <headnode-name> myScene.xml

The head node is responsible for all communication with your end, and it 
will in turn distribute work to the other cluster nodes. 

To get an overview of the currently running instances, run

$ ./cluster.py status

To shut all cluster nodes down (don't forget this step -- the running
instances will cost you money even when you are not using them), you can 
run

$ ./cluster.py terminateAll

In some cases (e.g. volume rendering with huge data files), it is impractical
to use Mitsuba in this way, because it will have to stream all scene data to the
cluster every time the scene is rendered. For this reason, there is a
'syncData' command, which simultaneously downloads one or more data files to
all cluster nodes (using Amazon's S3 service). You will have to install s3cmd
on your machine and do an initial connec to to S3, which will create a '.s3cfg'
configuration file in your home directory. Copy that file to the same
directory as the 'cluster.py' script. Afterwards, you should be able to run
e.g.

$ ./cluster.py syncData myHugeVolumeDataFile.vol

and the file will be available on all nodes shortly afterwards. You will
need to use special Mitsuba plugins, which can be made aware of local data 
files instead of streaming them (e.g. 'heterogeneous' or
'heterogeneous-flake' for volume rendering).
