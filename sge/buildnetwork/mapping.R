

## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

## read in all datasets
fin <- list()
k   <- 1
for( i in 1:4 ){
    PATHS[2] = sprintf("%s/%s",PATHS[1],datasets[i,1])
    st1      = sprintf("%s/%s.csv.gz",PATHS[2],datasets[i,2])
    if( file.exists(st1) && file.info(st1)$size!=0 ){
        fin[[k]] = data.table::fread(st1)
        names(fin)[k] = datasets[i,1]
        k=k+1
    }
}



# number of unique patient IDs across all datasets
ids <- c()
for( i in 1:length(fin) ){
    indx = which(colnames(fin[[i]])==IDS[ver,1])
    temp = as.character(unlist(fin[[i]][,.SD,.SDcols=indx[1]]))
    temp = temp[!duplicated(temp)]
    ids  = c(ids,temp)
}

ids = unique(ids)
write.table(as.data.frame(ids),sprintf("%s/%s.csv",PATHS[3],files[8]),sep="\t",row.names=F,col.names=F,quote=F)
#----


# build mapping list for each dataset
N = length(ids)

for( i in 1:length(fin) ){

    map=matrix(NA,nrow=N,ncol=4)
    colnames(map) <- c("Global patient Matrix indx","Global patient ID","Dataset patient ID","Dataset patient Matrix indx")
    map[,1] = seq_along(1:N)
    map[,2] = ids

    
    indx         = which(colnames(fin[[i]])==IDS[ver,1])    
    fids         = as.character(unlist(fin[[i]][,.SD,.SDcols=indx[1]]))
    findx        = match(fids,map[,2])
    map[findx,3] = fids
    map[findx,4] = seq_along(1:length(findx))

    map          = map[!is.na(map[,4]),]
    map          = map[order(as.numeric(map[,4])),]
    
    PATHS[2] = sprintf("%s/%s",PATHS[1],datasets[i,1])
    write.table(as.data.frame(map),sprintf("%s/%s_%s.csv",PATHS[2],datasets[i,2],files[7]),sep="\t",row.names=F,col.names=T,quote=F)

}
#----
