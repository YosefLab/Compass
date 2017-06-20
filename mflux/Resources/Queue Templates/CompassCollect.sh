#!/bin/sh

DATA=$1
MODEL=$2
MEDIA=$3
TEMPDIR=$4

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

compass --data $DATA --model $MODEL --media $MEDIA --temp-dir $TEMPDIR --output-dir . --collect
