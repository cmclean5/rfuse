#!/bin/sh

echo "Running on Eddie..."

WORKINGDIR=/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK/buildnetwork/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute_super.sh

name="featureSel_20"

chmod +x $EXE

START=1
END=20

for i in `seq $START $END`
do
    qsub -N $name -l h_rt=12:00:00 -l h_vmem=8G $EXE
done

#  done

echo "$0 done!"

exit 0
    
