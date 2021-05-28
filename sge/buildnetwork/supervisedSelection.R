## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

## load variable selection functions
source("loadVariableReductionFunctions.R")

source('loadMI.R')

## supervised learning libraries
#library(kernlab)      # SVM methodology
#library(e1071)        # SVM methodology
library(RRF)          # regularized random forest algorithm, based on the randomForest R package

args       <- commandArgs(TRUE);
ITS        <- as.numeric(args[1])
SEED       <- as.numeric(args[2])
ASDsubset  <- as.numeric(args[3])
gamma      <- as.numeric(args[4])

if( gamma < 0 ){ gamma = 0; } 
if( gamma > 1 ){ gamma = 1; }

cat("ITS         : ", ITS, "\n")
cat("SEED        : ", SEED, "\n")
cat("ASDsubset   : ", ASDsubset, "\n")
cat("gamma       : ", gamma, "\n")

## set random number seed
setSeed(SEED);


##PATHS[1]=DIRS[2]
PATHS[2]=sprintf("%s/%s",DIRS[1],"buildnetwork")

## set MI data directory
SubConDir = "thres9"
SubASDDir = "thres7"

ExtCon   = "fuse_Pc_control"
ExtASD   = "fuse_Pc"

Ext2Con    = "control"
Ext2ASD    = "asd"

DSConDir <- sprintf("%s/datasets/%s/%s",DIRS[1],version[ver],"control")
DSASDDir <- sprintf("%s/datasets/%s",DIRS[1],version[ver])

MIConDir <- sprintf("%s/MATRIX/%s",PATHS[2],SubConDir)
MIASDDir <- sprintf("%s/MATRIX/%s",PATHS[2],SubASDDir)

##controlDir <- sprintf("%s/RANDOM/%s/",PATHS[2],SubConDir)
##ASDDir     <- sprintf("%s/RANDOM/%s/",PATHS[2],SubASDDir)


## set files specific to addPatientAnnnotation.R script
files <- matrix(NA,ncol=2,nrow=5)
files[1,1] = "ind";      files[1,2] = "individuals"
files[2,1] = "bms";      files[2,2] = "basic_medical_screening"
files[3,1] = "scq";      files[3,2] = "scq"
files[4,1] = "rbsr";     files[4,2] = "rbsr"
files[5,1] = "dcdq";     files[5,2] = "dcdq"
#----

con <- Datastudy( DIR=PATHS[2], DSDir=DSConDir, MIDir=MIConDir, SubDir=SubConDir, Ext=ExtCon, ASDsubset=0 )

asd <- Datastudy( DIR=PATHS[2], DSDir=DSASDDir, MIDir=MIASDDir, SubDir=SubASDDir, Ext=ExtASD, ASDsubset=ASDsubset )

z1 = con$ResC
z2 = asd$ResC

## deal with NA
z1 = dealWithNA(MAT=z1,NAval=-1)
z2 = dealWithNA(MAT=z2,NAval=-1)

rownames(z1) = rep("con",length(z1[,1]))
rownames(z2) = rep("asd",length(z2[,1]))

data = rbind(z1,z2)
Y    = c(rownames(z1),rownames(z2))   

indx_train   <- list()
trainData    <- list()
testData     <- list()
trainY       <- list()
testY        <- list()
RF           <- list()
GRF          <- list()
GRF_RF       <- list()
predRF       <- list()
fscoreRF     <- list()
predGRF      <- list()
fscoreGRF    <- list()
predGRF_RF   <- list()
fscoreGRF_RF <- list()

for( p in 1:ITS ){

     cat("run ", p,"...")

     indx_train[[p]] <- sample.int(nrow(data), size = 0.75*nrow(data))

     trainData[[p]]  <- data[indx_train[[p]],];  trainY[[p]] <- Y[indx_train[[p]]];
     testData[[p]]   <- data[-indx_train[[p]],]; testY[[p]]  <- Y[-indx_train[[p]]];

     ## build an ordinary random forest
     RF[[p]] <- RRF( trainData[[p]], flagReg=0, as.factor(trainY[[p]]))##, keep.forest=TRUE )

     imp   <- RF[[p]]$importance[,"MeanDecreaseGini"]
     impRF <- imp/max(imp) # normalization

     ## build a guided random forest (GRF) with gamma = 1 (weighting applied to feature's importance score,
     ## gamma = 1 implies maximum penalty on less important features, leading to small feature set).
     ## note the difference between GRF and RF is that ’coefReg’ is related to impRF in GRF, while
     ## it is constant for all variables in RF.
     ## gamma   <- 0.5

     coefReg <- (1-gamma) +gamma*impRF

     ## build guided random forest
     GRF[[p]]     <- RRF(trainData[[p]],as.factor(trainY[[p]]),flagReg=0,coefReg=coefReg)##, keep.forest=TRUE)

     ## build random forest using guided random forest's (GRF) feature subset
     GRF_RF[[p]]<-RRF(trainData[[p]][,GRF[[p]]$feaSet],flagReg=0,as.factor(trainY[[p]]))##,keep.forest=TRUE)


     ## predict on test data using RF
     predRF[[p]]           <- predict(RF[[p]],testData[[p]])
     fscoreRF[[p]]         <- fscore(cbind(testY[[p]],as.vector(predRF[[p]])))
     cat("*** RF fscore ****")
     print(fscoreRF[[p]])
     cat("*******")

     ## predict on test data using GRF
     predGRF[[p]]   <- predict(GRF[[p]],testData[[p]])
     fscoreGRF[[p]] <- fscore(cbind(testY[[p]],as.vector(predGRF[[p]])))
     cat("*** GRF fscore ****")
     print(fscoreGRF[[p]])
     cat("*******")

     ## predict on test data using GRF_RF
     predGRF_RF[[p]]   <- predict(GRF_RF[[p]],testData[[p]][,GRF[[p]]$feaSet])
     fscoreGRF_RF[[p]] <- fscore(cbind(testY[[p]],as.vector(predGRF[[p]])))
     cat("*** GRF_RF fscore ****")
     print(fscoreGRF_RF[[p]])
     cat("*******")

     cat("done.\n")
}

save(indx_train,trainData,testData,trainY,testY,RF,GRF,GRF_RF,SEED,gamma,ITS,
     predRF,fscoreRF,predGRF,fscoreGRF,predGRF_RF,fscoreGRF_RF, file = "results.RData"); 

##load(“my_data.RData”)

## NOTES:
## with gamma == 1,   only scq_final_score selected.
## with gamma == 0.8, scq_final_score, scq_q39_imaginative_games
## with gamma == 0.2, 123 features selected
## with gamma == 0.5

##-----Example 1 -----
#library(RRF);set.seed(1)
#     
##only the first feature and last feature are truly useful
#X <- matrix(runif(50*50), ncol=50)
#class <- (X[,1])^2 + (X[,50])^2  
#class[class>median(class)] <- 1;
#class[class<=median(class)] <- 0
#     
##ordinary random forest. 
#rf <- RRF(X,as.factor(class), flagReg = 0)
#impRF <- rf$importance
#impRF <- impRF[,"MeanDecreaseGini"]
#rf$feaSet
#     
##regularized random forest
#rrf <- RRF(X,as.factor(class), flagReg = 1)
#rrf$feaSet
#     
##guided regularized random forest
#imp <- impRF/(max(impRF))#normalize the importance score
#gamma <- 0.5
#coefReg <- (1-gamma)+gamma*imp #weighted average
#grrf <- RRF(X,as.factor(class),coefReg=coefReg, flagReg=1)
#grrf$feaSet
#     
##guided random forest
#gamma <- 1
#coefReg <- (1-gamma)+gamma*imp 
#grf <- RRF(X,as.factor(class),coefReg=coefReg, flagReg=0)
#grf$feaSet
 
