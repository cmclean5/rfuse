## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

args       <- commandArgs(TRUE);
CORES      <- as.numeric(args[1])
LOCALPATH  <- as.character(args[2])
RUNS       <- as.numeric(args[3])
NR         <- as.numeric(55905)

if( CORES != nCORES ){
    nCORES = CORES
}

## probably need to change this, but fuseDataset will always perform 1 run, 
## so first remove 1 from RUNS. 
RUNS = RUNS-1

if( RUNS <= 0 ) { RUNS = 0; }

if( RUNS >= 10 ){ RUNS = 9; }

cat("nCORES   : ", nCORES, "\n")
cat("LOCALPATH: ", LOCALPATH,"\n") 	
cat("RUNS     : ", RUNS, "\n")
cat("NR       : ", NR,"\n") 	

run <- vector(length=2)
run[1] = 1 # copy and load global matrices for each dataset
run[2] = 1 # fuse global matrices from each dataset

#----------------------------------
# done
#----------------------------------

#----------------------------------
# initialise new fusion model (fuse.cpp)
#----------------------------------

cat("#----------------------------------\n")
cat("#  create fusion model             \n")
cat("#----------------------------------\n")     

#create model first
cat("creating model..."); rfuse::createModel(); cat("done.\n");
#------

cat("#----------------------------------\n")
cat("# done                             \n")
cat("#----------------------------------\n")

#----------------------------------
# done
#----------------------------------


if( run[1] ){

    #----------------------------------
    # copy and load matrices for each dataset
    #----------------------------------

    cat("#----------------------------------\n")
    cat("# copy and load matrices for each dataset\n")
    cat("#----------------------------------\n")       

    cat("#----------------------------------\n")
    cat("# copy bms dataset\n")
    cat("#----------------------------------\n") 

    #bms
    g=1
    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
    PATHS[5]=sprintf("%s/MATRIX/%s",PATHS[2],ashort[thres])

    copyDataset( DATASET=datasets[g,1], PATH=PATHS[5] )
    changeBK   ( DATASET=datasets[g,1], LOCALPATH=LOCALPATH )

    Pbms  = readRDS(sprintf("%s_Pr.rds",datasets[g,1]))    

    ##record the global matrix size
    NR    = as.numeric(Pbms$nrow)
    rm(Pbms); gc();
    
    loadDataset( DATASET=datasets[g,1], SET=c(1), nCORES=nCORES );

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")

    cat("#----------------------------------\n")
    cat("# copy scq dataset\n")
    cat("#----------------------------------\n") 

    #scq
    g=2
    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
    PATHS[5]=sprintf("%s/MATRIX/%s",PATHS[2],ashort[thres])

    copyDataset( DATASET=datasets[g,1], PATH=PATHS[5] )
    changeBK   ( DATASET=datasets[g,1], LOCALPATH=LOCALPATH )
    loadDataset( DATASET=datasets[g,1], SET=c(2), nCORES=nCORES );

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")

    if( buildASD == TRUE ){

    cat("#----------------------------------\n")
    cat("# copy rbsr dataset\n")
    cat("#----------------------------------\n") 

    #rbsr
    g=3
    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
    PATHS[5]=sprintf("%s/MATRIX/%s",PATHS[2],ashort[thres])

    copyDataset( DATASET=datasets[g,1], PATH=PATHS[5] )
    changeBK   ( DATASET=datasets[g,1], LOCALPATH=LOCALPATH )
    loadDataset( DATASET=datasets[g,1], SET=c(3), nCORES=nCORES );

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")

    cat("#----------------------------------\n")
    cat("# copy dcdq dataset\n")
    cat("#----------------------------------\n") 

    #dcdq
    g=4
    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
    PATHS[5]=sprintf("%s/MATRIX/%s",PATHS[2],ashort[thres])

    copyDataset( DATASET=datasets[g,1], PATH=PATHS[5] )
    changeBK   ( DATASET=datasets[g,1], LOCALPATH=LOCALPATH )
    loadDataset( DATASET=datasets[g,1], SET=c(4), nCORES=nCORES );

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")
    
   }

    
}#run1


if( run[2] ){

    #----------------------------------
    # fuse global matrices for each dataset
    #----------------------------------

    cat("#----------------------------------\n")
    cat("# fuse global matrices from each dataset\n")
    cat("#----------------------------------\n")
   
    ## create the global fused probability matrix on disk and link to it.
    Pc  =  FBM(NR,NR,backingfile=file.path(LOCALPATH,fuseName[1]),FBMtype[FBMt])
    #----                     
    
    cat("Start fusion algorithm...\n");

    if( buildASD == TRUE ){
        cat("build ASD network, \n");

	nTHREADS=nCORES

        res = rfuse::fuseASDsets(Pc, nCORES, nTHREADS, RUNS );
        rfuse::normalisePc( Pc, nCORES )        

    } else {
        cat("build control network, \n");     
	
        nTHREADS=nCORES

        res = rfuse::fuseCONTROLsets(Pc, nCORES, nTHREADS, RUNS );
	rfuse::normalisePc( Pc, nCORES )           
    
    }      
   
    cat("done.\n")
  

    ## write results table
    write.table(res$RESULT,sprintf("%s/%s.csv",LOCALPATH,fuseName[1]),row.names=,col.names=T,quote=F)
    
    ## save the global fused probability matrix
    Pc$save();
          
    rm(Pc); gc();

    ## new output path for Pc
    PATHS[3]=sprintf("%s/%s",DIRS[1],"buildnetwork")
    PATHS[5]=sprintf("%s/MATRIX/%s",PATHS[3],ashort[thres])

    ## copy Pc to output dir
    copy2OutDir( DATASET=fuseName[1], LOCALPATH=LOCALPATH, OUTPATH=PATHS[5] );       

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")
    
    
}#run2

#----------------------------------
# delete fusion model (fuse.cpp)
#----------------------------------

cat("#----------------------------------\n")
cat("#  delete fusion model             \n")
cat("#----------------------------------\n")     

#delete model
cat("deleting model..."); rfuse::deleteModel(); cat("done.\n");

cat("#----------------------------------\n")
cat("# done                             \n")
cat("#----------------------------------\n")


#----------------------------------
# done
#----------------------------------
