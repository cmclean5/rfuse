#!/bin/sh

echo "Running on Eddie..."

WORKINGDIR=/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK/buildnetwork/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

#nCORES=10
nCORES=12

nRUNS=0
##nRUNS=1
#nRUNS=5
##nRUNS=10

#h_vmem==25G,30G

name="fuse_20"

chmod +x $EXE

#for i in `seq $START $END`
#    do

##nRUNS=0
qsub -N $name -l h_rt=01:00:00 -pe sharedmem $nCORES -l h_vmem=30G -v nCORES=$nCORES -v nRUNS=$nRUNS $EXE

##nRUNS=1
##qsub -N $name -l h_rt=04:00:00 -pe sharedmem $nCORES -l h_vmem=30G -v nCORES=$nCORES -v nRUNS=$nRUNS $EXE

##nRUNS=5
#qsub -N $name -l h_rt=24:00:00 -pe sharedmem $nCORES -l h_vmem=30G -v nCORES=$nCORES -v nRUNS=$nRUNS $EXE

##nRUNS=10
##qsub -N $name -l h_rt=24:00:00 -pe sharedmem $nCORES -l h_vmem=30G -v nCORES=$nCORES -v nRUNS=$nRUNS $EXE

#  done

echo "$0 done!"

exit 0
    
