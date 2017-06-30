#!/bin/sh

CONFIG=$1

workdir=$PBS_O_INITDIR

mkdir -p $workdir
exec 2> $workdir/err.log > $workdir/out.log

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
compass --data '' --collect --config-file $CONFIG
