## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

Nmat = 500

ii   = seq(1,(Nmat+1),1)
II   = cbind(ii,(ii+1))
II   = cbind(II,rep(0,length(II[,1])),rep(NA,length(II[,1])) )


#files    <- vector(length=2)
#files[1] <- "indices.csv"
#files[2] <- "result.csv.gz"

## reset paths for missingIndices.R script
PATHS[2] = sprintf("%s/%s",PATHS[1],datasets[choice,1])
PATHS[3]=sprintf("%s/EDDIE/RESULTS",PATHS[2])
##---

subdirs  = list.files(path=sprintf("%s",PATHS[3]));
nstudies = length(subdirs);


for( f in 1:nstudies ){

    st1 = sprintf("%s/%s/%s",PATHS[3],subdirs[f],files[1]);
    st2 = sprintf("%s/%s/%s",PATHS[3],subdirs[f],files[2]);	
    if( file.exists(st1) && file.info(st1)$size!=0 &&
        file.exists(st2) && file.info(st2)$size!=0 ){
     tb = read.table(st1,header=F,sep="\t");     
     indx = which(II[,1]==tb[[1]][1] & II[,2]==tb[[1]][2])
     if( length(indx) != 0 ){ II[indx[1],3] = 1; II[indx[1],4] = subdirs[f]; }
    }

     rm(tb)
}


reruns = II[as.numeric(II[,3]) == 0,2]

## write 'summary.csv'
write.table(II,sprintf("%s/%s.csv",PATHS[2],files[3]),sep="\t",col.names=F,row.names=F,quote=F)

## write 'reruns.csv'
write.table(reruns,sprintf("%s/%s.csv",PATHS[2],files[4]),sep="\t",col.names=F,row.names=F,quote=F)

