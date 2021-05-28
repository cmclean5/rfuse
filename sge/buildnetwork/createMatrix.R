

## run setParamters.R script
source("setParameters.R")

## run setUp.R script
source("setUp.R")

## set specific paths in createMatrix.R
PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[choice,1])
PATHS[5]=sprintf("%s/MATRIX",PATHS[2])

## compile c++ code -- 'bigstatsUtils2.cpp' file.
sourceCpp(sprintf("%s/%s",PATHS[4],files[5]))

## read in dataset filtered file
fin   = data.table::fread(sprintf("%s/%s.csv.gz",PATHS[2],datasets[choice,2]))
## read in dataset 'summary.csv' file
zz    = read.delim(sprintf("%s/%s.csv",PATHS[2],files[3]),sep="\t",header=F)

NR = as.numeric(length(rownames(fin)))
NC = as.numeric(length(colnames(fin)))

## create a big matrix which is back on disk

## create temp dir: "/tmp/RtmpXXXX/"
#td <- tempdir()
td <- PATHS[5]
print(td)

## clean space, if matrix already exist on disk. 
st1 = sprintf("%s/%s.bk",td,datasets[choice,1])
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }
st1 = sprintf("%s/%s.rds",td,datasets[choice,1])
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }

## create the matrix on disk and link to it.
mat = FBM(NR,NR,backingfile=file.path(td,datasets[choice,1]), FBMtype[FBMt] )
    
# Start the clock!
ptm <- proc.time()

for( i in 1:length(zz[,1]) ){

    ## read in EDDIE 'result.csv.gz' files
    st1 = sprintf("%s/EDDIE/RESULTS/%s/%s",PATHS[2],as.character(zz[i,4]),files[2]);

    if( file.exists(st1) && file.info(st1)$size!=0 ){
        sim      = data.table::fread(st1)
        nn       = as.numeric(length(unlist(sim[,1])))
        if( nn > 0 ){

            if( i %% 10 == 0 ){ cat("output: ", round((i/length(zz[,1]))*100, 3), "% done.\n") }
            
            indxi = as.numeric(unlist(sim[,1]))
            indxj = as.numeric(unlist(sim[,2]))
            val   = as.numeric(unlist(sim[,13]))      

            #ele   = (indxj-1)*NR+(indxi-1)+1
            #rele   = (indxi-1)*NR+(indxj-1)+1
            #mat[ele]  = val
            #mat[rele] = val    
            
            addEntries(mat,indxi,indxj,val)
                

        }
    }
}
    
# Stop the clock
pet <- proc.time() - ptm
cat(sprintf("time = %.3f \n", pet[[1]]))

## create rds files
mat$save()
    


