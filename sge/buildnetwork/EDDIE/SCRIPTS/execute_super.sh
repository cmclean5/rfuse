#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK/buildnetwork

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

# load scripts
cp -r $EXECDIR/setUp.R .
cp -r $EXECDIR/setParameters.R .
cp -r $EXECDIR/loadFunctions.R .
cp -r $EXECDIR/loadVariableReductionFunctions.R .
cp -r $EXECDIR/loadMI.R .
cp -r $EXECDIR/supervisedSelection.R .

# initialise environment module
. /etc/profile.d/modules.sh

#load module R
module load igmm/apps/R/3.6.3 #for rfuse library

## tell the Unix shell to use the maximum amount of stack memory available
ulimit -s unlimited

#print memory limit
echo "MEMORY LIMIT... "
ulimit -v
echo "...done."
#----------

## find candiate edges in SPARK datasets
time Rscript supervisedSelection.R 10 $JOB_ID 1 1

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.RData $OUTDIR
cp -v *.RDS   $OUTDIR
cp -v *.csv   $OUTDIR

echo "$0 done!"

exit 0
