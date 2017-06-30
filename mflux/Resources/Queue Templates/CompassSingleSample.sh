#!/bin/sh

CONFIG=$1

workdir=$PBS_O_INITDIR/sample$PBS_ARRAYID

mkdir -p $workdir
exec 2> $workdir/err.log > $workdir/out.log

## Remove the old log files
rm $PBS_O_HOME/$PBS_JOBNAME.o*
rm $PBS_O_HOME/$PBS_JOBNAME.e*

### Switch to the working directory
echo Working directory is $workdir
cd $workdir

### Run some informational commands.
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`

# --data is taken from the config file
compass --data '' --single-sample $PBS_ARRAYID --config-file $CONFIG
