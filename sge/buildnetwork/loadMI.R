MIstudy <- function( DIR=NULL, DSDir=NULL, MIDir=NULL, SubDir=NULL, Ext=NULL, Ext2=NULL, ASDsubset=NULL, Ychoice=NULL ){

    ## DIR       == PATHS[2]
    ## DSDir     == PATHS[3]
    ## MIDir     == MIdir
    ## SubDir    == ashort[thres]
    ## Ext       == ggName[1]
    ## Ext2      == extName[1]
    ## ASDsubset == ASDsubset
    ## Ychoice   == Ychoice
    
    ## read in datasets
    fin  <- list()
    temp <- c()

    Nfiles = dim(files)[1]

    
    ## select features in ASD which match those in control
    if( ASDsubset ){ Nfiles = 2; }

    for( i in 1:Nfiles ){

        if( dir.exists(DSDir) ){    
            st1 = sprintf("%s/%s/%s.csv.gz",DSDir,datasets[i,1],datasets[i,2])
            if( file.exists(st1) && file.info(st1)$size!=0 ){
                fin[[i]] = data.table::fread(st1,sep=SEP[1],header=T)
                setkeyv(fin[[i]],IDS[1])
                names(fin)[i] <- datasets[i,1]
                temp          <- c(temp,as.vector(unlist(fin[[i]][,1])))
            }
        }
        
    }


    vatts = read.delim(sprintf("%s/%s/%s/%s.csv",DIR,"VATTS",SubDir,Ext),sep="\t",header=T)

    ResC <- buildPhenotypeMatrix( VATTS=vatts, FIN=fin, ID="SPARKID" )

    
## define class labels ##
Y1   = vatts[,which(colnames(vatts)=="fcSc")]  # patient clustering at level 1
Y2   = vatts[,which(colnames(vatts)=="fcSc5")] # patient clustering at level 5

Cind    <- which(colnames(vatts)=="GENDER")
Classes <- levels(vatts[,Cind])
Y3   = match(vatts[,Cind],Classes) # Gender
 
Cind    <- which(colnames(vatts)=="DIAG")
Classes <- levels(vatts[,Cind])
Y4   = match(vatts[,Cind],Classes) # DIAG 

Cind    <- which(colnames(vatts)=="DIAG_AGE_MONTHS")
Classes <- levels(as.factor(vatts[,Cind]))
Y5   = match(vatts[,Cind],Classes) # DIAG AGE MONTHS

#we'll get many classes using month, i.e. class imbalance, so convert months to years.
Y5   = ceiling(Y5/12)
    
Cind    <- which(colnames(vatts)=="CONG")
Classes <- levels(vatts[,Cind])
Y6   = match(vatts[,Cind],Classes) # cognitive impairment diagnosis
 
Cind    <- which(colnames(vatts)=="LANG")
Classes <- levels(vatts[,Cind])
Y7   = match(vatts[,Cind],Classes) # LANGUAGE LEVEL

Cind    <- which(colnames(vatts)=="rGRADE")
Classes <- levels(vatts[,Cind])
Y8   = match(vatts[,Cind],Classes) # repeat GRADE

Cind    <- which(colnames(vatts)=="dADULT")
Classes <- levels(vatts[,Cind])
Y9   = match(vatts[,Cind],Classes) # dependent adult

Cind    <- which(colnames(vatts)=="mBIRTH")
Classes <- levels(vatts[,Cind])
Y10   = match(vatts[,Cind],Classes) # multiple birth

Cind    <- which(colnames(vatts)=="ENROL")
Classes <- levels(vatts[,Cind])
Y11   = match(vatts[,Cind],Classes) # multiple enrolled

Cind    <- which(colnames(vatts)=="RACE")
Classes <- levels(vatts[,Cind])
Y12   = match(vatts[,Cind],Classes) # multiple enrolled
    
X       <- printPer6Months(vatts, "", "BOWEL", PRINT=FALSE)
Classes <- levels(as.factor(X))
Y13   = match(X,Classes) # bowel trained

X       <- printPer6Months(vatts, "", "cPHRASES", PRINT=FALSE)
Classes <- levels(as.factor(X))
Y14   = match(X,Classes) # complex phrases

    
if( buildASD ){
    Y = cbind(Y1,Y3,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y4,Y5,Y6,Y7)
    colnames(Y) <- c("level1","GENDER","rGRADE","dADULT","mBIRTH","ENROL","RACE","BOWEL","cPHRASES","DIAG","DIAG_AGE_YEARS","CONG","LANG")
    Yweight <- rep(1,dim(Y)[2]) ##each class should have a weighting [0,1]
} else {
    Y = cbind(Y1,Y3,Y8,Y9,Y10,Y11,Y12,Y13,Y14)
    colnames(Y) <- c("level1","GENDER","rGRADE","dADULT","mBIRTH","ENROL","RACE","BOWEL","cPHRASES")
    Yweight <- rep(1,dim(Y)[2]) ##each class should have a weighting [0,1]
}
    
  
##----

## step 1: select subset of class labels
Ysub           = as.data.frame(Y[,Ychoice])
YsubWe         = Yweight[Ychoice]
colnames(Ysub) = colnames(Y)[Ychoice]
Nclasses       = length(colnames(Ysub))


## step 2: build feature MI given ALL class labels
su = symmetricUncertainty( MAT=ResC, GROUPS=Ysub )

## step 3: load feature MI given each class label
res       <- list()
NMIxixkCY <- list()
fNMI      <- list()

for( i in 1:Nclasses ){
    st1 = sprintf("%s/%s_%s.RDS",MIDir,colnames(Ysub)[i],Ext2)
    if( file.exists(st1) && file.info(st1)$size!=0 ){        
        res[[i]]      = readRDS(st1)
        names(res)[i] = colnames(Ysub)[i]
        
        NMIxixkCY[[i]]      = res[[i]]$NMI
        if( ASDsubset ){ NMIxixkCY[[i]] = extractSubMatrix( NMIxixkCY[[i]], colnames(ResC) ); }
        names(NMIxixkCY)[i] = names(res)[i]
    }
}
    
if( ASDsubset ){
    for( i in Nclasses ){
        res[[i]]$MI     = extractSubMatrix(res[[i]]$MI,     colnames(ResC) );
        res[[i]]$NMI    = extractSubMatrix(res[[i]]$NMI,    colnames(ResC) );
        res[[i]]$II     = extractSubMatrix(res[[i]]$II,     colnames(ResC) );
        res[[i]]$IXY    = extractSubMatrix(res[[i]]$IXY,    colnames(ResC) );
        res[[i]]$NMIXY  = extractSubMatrix(res[[i]]$NMIXY,  colnames(ResC) );
        res[[i]]$IXZ    = extractSubMatrix(res[[i]]$IXZ,    colnames(ResC) );
        res[[i]]$NMIXZ  = extractSubMatrix(res[[i]]$NMIXZ,  colnames(ResC) );
    }
}

## step 4: load feature MI
st1 = sprintf("%s/%s_%s.RDS",MIDir,"featureNMI",Ext2)
if( file.exists(st1) && file.info(st1)$size!=0 ){
    fNMI = readRDS(st1)
    if( ASDsubset ){ fNMI = extractSubMatrix(fNMI, colnames(ResC) ); }
}


## expected feature correlation condition on class
CNMIavg = matrix(0,nrow=dim(NMIxixkCY[[1]])[1],ncol=dim(NMIxixkCY[[1]])[2])
colnames(CNMIavg) = colnames(NMIxixkCY[[1]])
rownames(CNMIavg) = rownames(NMIxixkCY[[1]])
        
for( i in 1:Nclasses ){ CNMIavg = CNMIavg + YsubWe[i] * NMIxixkCY[[i]]; }

CNMIavg = CNMIavg/Nclasses; 
    
    return(list(ResC=ResC,Y=Y,su=su,res=res,NMIxixkCY=NMIxixkCY,CNMIavg=CNMIavg, fNMI=fNMI))
    
}

Datastudy <- function( DIR=NULL, DSDir=NULL, MIDir=NULL, SubDir=NULL, Ext=NULL, ASDsubset=NULL ){

    ## DIR       == PATHS[2]
    ## DSDir     == PATHS[3]
    ## MIDir     == MIdir
    ## SubDir    == ashort[thres]
    ## Ext       == ggName[1]
    ## ASDsubset == ASDsubset
    
    ## read in datasets
    fin  <- list()
    temp <- c()

    Nfiles = dim(files)[1]

    
    ## select features in ASD which match those in control
    if( ASDsubset ){ Nfiles = 2; }

    for( i in 1:Nfiles ){

        if( dir.exists(DSDir) ){    
            st1 = sprintf("%s/%s/%s.csv.gz",DSDir,datasets[i,1],datasets[i,2])
            if( file.exists(st1) && file.info(st1)$size!=0 ){
                fin[[i]] = data.table::fread(st1,sep=SEP[1],header=T)
                setkeyv(fin[[i]],IDS[1])
                names(fin)[i] <- datasets[i,1]
                temp          <- c(temp,as.vector(unlist(fin[[i]][,1])))
            }
        }
        
    }


    vatts = read.delim(sprintf("%s/%s/%s/%s.csv",DIR,"VATTS",SubDir,Ext),sep="\t",header=T)

    ResC <- buildPhenotypeMatrix( VATTS=vatts, FIN=fin, ID="SPARKID" )

    return(list(ResC=ResC))
    
}
