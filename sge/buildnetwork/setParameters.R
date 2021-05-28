#----------------------------------
# set parameter values
# 
# Note: building a Minimum Spanning Tree MST,
# see this link: https://rpubs.com/phat/mst
#----------------------------------


#--- reset
rm(list=ls())

## bigstatsr FBM data type
FBMt=1 # 1 == "double", 2 == "float"

## are we builing a ASD or control network
#buildASD=TRUE
buildASD=FALSE
#----

## set SPARK dataset version
#ver=1 # ==> 'c3'
ver=2 # ==> 'c4'
#----

## set the choice of dataset 
choice=1 # ==> 'bms'
#choice=2 # ==> 'scq'
#choice=3 # ==> 'rbsr'
#choice=4 # ==> 'dcdq'
#choice=5 # ==> 'joined'
#-----

#-------------------------------------
# set thres parameter in 'buildnetwork.R' script
#-------------------------------------
#
#-------------------------------------
# Notes: for building networks from SPARK version 'c3' datasets
#      : ADS network was built setting thres=5, i.e. p<=0.1,
#      : control network was built setting thres=1, i.e. p<=0.5
#      : initially, for joined ASD network set thres=7, i.e. p<=0.01)
#-------------------------------------
#
#-------------------------------------
# Notes: for building networks from SPARK version 'c4' datasets
#      : ADS network was built setting thres=6, i.e. p<=0.05,
#      : control network was built setting thres=1, i.e. p<=0.5
#      : initially, for joined ASD network set thres=5, i.e. p<=0.01)
#-------------------------------------
#
#
thres=9
#----

#----------------------------------
# done
#----------------------------------

#-------------------------------------
# Order for running scripts
# 1) missingIndices.R
# 2) createMatrix.R
# 3) uniqueMatrixes.R (check for any duplicate rows in matrix) 
# 3) buildnetwork.R (using option run[1]=1)
# 4) mapping.R      (using option run[1]=1)
# 5) mapping.R      (using option run[2]=1)
#-------------------------------------

