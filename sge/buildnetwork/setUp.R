#----------------------------------
# setup
#----------------------------------

## libraries
library(data.table)
library(Rcpp)
library(bigstatsr)
library(igraph)
library(R.utils)
library(ggplot2)
library(rfuse)

## SPARK datasets
datasets <- matrix(NA,ncol=2,nrow=5)
datasets[1,1] = "bms";    datasets[1,2] = "basic_medical_screening_filtered";
datasets[2,1] = "scq";    datasets[2,2] = "scq_filtered";
datasets[3,1] = "rbsr";   datasets[3,2] = "rbsr_filtered";
datasets[4,1] = "dcdq";   datasets[4,2] = "dcdq_filtered";
datasets[5,1] = "joined"; datasets[5,2] = "patient_joinedData";

## version of the SPARK datasets
version <- vector(length=2)
version[1]="c3"
version[2]="c4"

#---Directories
DIRS    <- vector(length=3)
DIRS[1] <- "/exports/csce/eddie/inf/groups/statbio/cmclean/SPARK"
DIRS[2] <- "/home/cmclean/WORK/DATA/SimonsData240120/SPARK_Collection_Version3"
DIRS[3] <- "/home/cmclean/WORK/DATA/SimonsData010720/SparkRapidDataRelease2020-07-01"

## files
files <- vector(length=8)
files[1]="indices.csv"
files[2]="result.csv.gz"
files[3]="summary"
files[4]="reruns"
files[5]="bigstatsUtils2.cpp"
files[6]="fuse.cpp"
files[7]="mapping"
files[8]="patientIDs_acrossDatasets"

## patient ids
IDS <- matrix(NA,ncol=2,nrow=2)
IDS[1,1] = "subject_sp_id"; IDS[2,1] = "subject_sp_id";
IDS[1,2] = "family_id";     IDS[2,2] = "family_sf_id";

## are we builing a ASD or control network
if( buildASD == TRUE ){

    PATHS <- vector(length=6)
    PATHS[1]=sprintf("%s/datasets/%s",DIRS[1],version[ver])
    PATHS[2]=""
    PATHS[3]=sprintf("%s/%s",DIRS[1],"buildnetwork")
    PATHS[4]=sprintf("%s/%s",PATHS[3],"rcpp")
    PATHS[5]=""
    PATHS[6]=""   

    SCQ <- matrix(NA,ncol=2,nrow=2)
    SCQ[1,1] = 10; SCQ[1,2] = 10;
    SCQ[2,1] = 30; SCQ[2,2] = 30;

    RBSR <- matrix(NA,ncol=2,nrow=2)
    RBSR[1,1] = 10;  RBSR[1,2] = 10;
    RBSR[2,1] = 100; RBSR[2,2] = 100;
    
    DCDQ <- matrix(NA,ncol=2,nrow=2)
    DCDQ[1,1] = 20; DCDQ[1,2] = 20;
    DCDQ[2,1] = 60; DCDQ[2,2] = 60;

    fuseName    <- vector(length=1)
    fuseName[1] <- "fuse_Pc" 
    
    candEd      <- vector(length=1)
    candEd[1]   <- "CandEd"   

    ggName      <- vector(length=1)
    ggName[1]   <- fuseName[1];#"fuse_Pc"
     
    
} else { 
    
    PATHS <- vector(length=6)
    PATHS[1]=sprintf("%s/datasets/%s/control",DIRS[1],version[ver])
    PATHS[2]=""
    PATHS[3]=sprintf("%s/%s",DIRS[1],"buildnetwork")
    PATHS[4]=sprintf("%s/%s",PATHS[3],"rcpp")
    PATHS[5]=""
    PATHS[6]=""
    
    files[6]="fuse_control.cpp"
    files[8]="patientIDs_acrossDatasets_control"

    SCQ <- matrix(NA,ncol=2,nrow=2)
    SCQ[1,1] = 2;  SCQ[1,2] = 2;
    SCQ[2,1] = 30; SCQ[2,2] = 30;
      
    RBSR <- matrix(NA,ncol=2,nrow=2)
    RBSR[1,1] = 10;  RBSR[1,2] = 10;
    RBSR[2,1] = 100; RBSR[2,2] = 100;

    DCDQ <- matrix(NA,ncol=2,nrow=2)
    DCDQ[1,1] = 20; DCDQ[1,2] = 20;
    DCDQ[2,1] = 60; DCDQ[2,2] = 60;

    fuseName    <- vector(length=1)
    fuseName[1] <- "fuse_Pc_control" 
    
    candEd      <- vector(length=1)
    candEd[1]   <- "CandEd_control"

    ggName      <- vector(length=1)
    ggName[1]   <- fuseName[1];#"fuse_Pc_control"
   
}


## pvalue cut-off for building network in 'buildnetwork.R' script
alpha <- vector(length=10)
alpha[1]=0.5
alpha[2]=0.4
alpha[3]=0.3
alpha[4]=0.2
alpha[5]=0.1
alpha[6]=0.05
alpha[7]=0.01
alpha[8]=0.001
alpha[9]=0.15
alpha[10]=0.15

aname <- vector(length=10)
aname[1] = "0.5"
aname[2] = "0.4"
aname[3] = "0.3"
aname[4] = "0.2"
aname[5] = "0.1"
aname[6] = "0.05"
aname[7] = "0.01"
aname[8] = "0.001"
aname[9] = "0.15"
aname[10] = "0.15"

ashort <- vector(length=10)
ashort[1] = "thres1"
ashort[2] = "thres2"
ashort[3] = "thres3"
ashort[4] = "thres4"
ashort[5] = "thres5"
ashort[6] = "thres6"
ashort[7] = "thres7"
ashort[8] = "thres8"
ashort[9] = "thres9"
ashort[10] = "test"
#----

## fused matrices
#fuseName = vector(length=4)
#fuseName[1] = "fuse_Pc"         # with 0.1 threshold applied to datasets 
#fuseName[2] = "fuse_Pc_Sr2"
#fuseName[3] = "fuse_0.5_Pc"     # with 0.5 threshold applied to datasets 
#fuseName[4] = "fuse_0.5_Pc_Sr2"
#----

## SPARK class variables
SPclasses <- matrix(NA,ncol=2,nrow=15)
SPclasses[1,1]  = "sex";                        SPclasses[1,2]  = "sex";
SPclasses[2,1]  = "age_at_registration_months"; SPclasses[2,2]  = "age_at_registration_months";
SPclasses[3,1]  = "age_at_registration_years";  SPclasses[3,2]  = "age_at_registration_years";
SPclasses[4,1]  = "diagnosis_age";              SPclasses[4,2]  = "diagnosis_age";
SPclasses[5,1]  = "diagnosis";                  SPclasses[5,2]  = "diagnosis";
SPclasses[6,1]  = "cognitive_impairment";       SPclasses[6,2]  = "cognitive_impairment";
SPclasses[7,1]  = "language_level";             SPclasses[7,2]  = "language_level";
SPclasses[8,1]  = "gen_test_ep";                SPclasses[8,2]  = "gen_test_ep";
SPclasses[9,1]  = "gen_test_frax";              SPclasses[9,2]  = "gen_test_frax";
SPclasses[10,1] = "gen_test_id";                SPclasses[10,2] = "gen_test_id";
SPclasses[11,1] = "gen_test_mecp2";             SPclasses[11,2] = "gen_test_mecp2";
SPclasses[12,1] = "gen_test_nf1";               SPclasses[12,2] = "gen_test_nf1";
SPclasses[13,1] = "gen_test_noonan";            SPclasses[13,2] = "gen_test_noonan";
SPclasses[14,1] = "gen_test_pten";              SPclasses[14,2] = "gen_test_pten";
SPclasses[15,1] = "gen_test_tsc";               SPclasses[15,2] = "gen_test_tsc";
#----

## genetic testing
gtest    <- matrix(NA,ncol=2,nrow=8)
gtest[1,1] = "TEST_EPi";    gtest[1,2] = SPclasses[8,ver];#"gen_test_ep"
gtest[2,1] = "TEST_FRAX";   gtest[2,2] = SPclasses[9,ver];#"gen_test_frax"
gtest[3,1] = "TEST_ID";     gtest[3,2] = SPclasses[10,ver];#"gen_test_id"
gtest[4,1] = "TEST_MECP2";  gtest[4,2] = SPclasses[11,ver];#"gen_test_mecp2"
gtest[5,1] = "TEST_NF1";    gtest[5,2] = SPclasses[12,ver];#"gen_test_nf1"
gtest[6,1] = "TEST_NOONAN"; gtest[6,2] = SPclasses[13,ver];#"gen_test_noonan"
gtest[7,1] = "TEST_PTEN";   gtest[7,2] = SPclasses[14,ver];#"gen_test_pten"
gtest[8,1] = "TEST_TSC";    gtest[8,2] = SPclasses[15,ver];#"gen_test_tsc"
#----

##---clustering algorithms
alg    <- vector(length=2)
alg[1] <- "louvain"
alg[2] <- "fc"

##--- clustering algorithm weighted
algWe <- sprintf("%sWe",alg)

##--- clustering algorithm scaled
algSc <- sprintf("%sSc",alg)

##--- edge weighting schemes
we       <- vector(length=4)
we[1]    <- "NULL"
we[2]    <- "weight"
we[3]    <- "znorm"
we[4]    <- "scaled"


## fuse.cpp compile options
cxxflag.opts = "-I\"/usr/local/include/gperftools\" -I\"/usr/local/include/igraph\" -I\"/usr/local/lib/R/site-library/bigstatsr/include\" -DLEVEL1_DCACHE_LINESIZE=`getconf LEVEL1_DCACHE_LINESIZE` -O3 -march=native -mavx -fopenmp -ffast-math -DARMA_DONT_USE_WRAPPER"
#cxxflag.opts = "-I\"/usr/local/include/gperftools/\" -I\"/usr/local/include/igraph/\" -DLEVEL1_DCACHE_LINESIZE=`getconf LEVEL1_DCACHE_LINESIZE` -O3 -fopenmp -DARMA_DONT_USE_WRAPPER -lopenblas"

#lib.opts = "-fopenmp -L\"/usr/lib/x86_64-linux-gnu/\" -L\"/usr/local/lib/\" -ltcmalloc_minimal -ligraph"
lib.opts = "-L\"/usr/local/lib/\" -fopenmp -ltcmalloc_minimal -ligraph -lopenblas"

## bigstatsr data type
FBMtype <- vector(length=2)
FBMtype[1] <- "double"
FBMtype[2] <- "float"

## number of cores to use OpenMP
nCORES=c(4)

## for numerical stability   
SMALL = 1e-100

# define constants
EPSILON = .Machine$double.eps
#----

## plot dimensions
WIDTH=480
HEIGHT=480
#----

## delimers 
SEP <- vector(length=2)
SEP[1]="\t"
SEP[2]=","
#----

#----------------------------------
# load utility functions
#----------------------------------

source("loadFunctions.R")

#----------------------------------
# done
#----------------------------------
