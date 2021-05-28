#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK/buildnetwork
#SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/SBM/EDDIE/SCRIPTS
#DATADIR=/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK/datasets/c4/bms/MATRIX/thres6

CORES=$nCORES
RUNS=$nRUNS

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

# load scripts
cp -r $EXECDIR/setUp.R .
cp -r $EXECDIR/setParameters.R .
cp -r $EXECDIR/loadFunctions.R .
cp -r $EXECDIR/fuseMatrices.R .
cp -r $EXECDIR/getEdges.R .

# initialise environment module
. /etc/profile.d/modules.sh

#openmp environment
export OMP_DISPLAY_ENV=false
#export OMP_STACKSIZE=512M
export OMP_STACKSIZE=10G
export OMP_WAIT_POLICY=active 
export OMP_DYNAMIC=false 
export OMP_PROC_BIND=true #try to prevent the operating system from migrating a thread . This prevents many scaling problems.
export OMP_PLACES=threads

echo "OPENMP ENVIRONMENT... "
echo "OMP_DISPLAY_ENV: " $OMP_DISPLAY_ENV
echo "OMP_STACKSIZE:   " $OMP_STACKSIZE
echo "OMP_WAIT_POLICY: " $OMP_WAIT_POLICY
echo "OMP_DYNAMIC:     " $OMP_DYNAMIC
echo "OMP_PROC_BIND:   " $OMP_PROC_BIND
echo "OMP_PLACES:      " $OMP_PLACES
echo "...done."
#--------

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
time Rscript getEdges.R $CORES $TMPDIR $RUNS

## main fusion code, fuse SPARK datasets
#time Rscript fuseMatrices.R $CORES $TMPDIR $RUNS

#OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

#if [ ! -d $OUTDIR ]; then
#    mkdir $OUTDIR
#fi

#cp -v *.RDS $OUTDIR
#cp -v *.csv $OUTDIR

echo "$0 done!"

exit 0
