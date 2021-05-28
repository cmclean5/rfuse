## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

#----------------------------------
# done
#----------------------------------


#----------------------------------
# check any duplicate rows in matrices for each dataset
#----------------------------------

cat("#----------------------------------\n")
cat("# check for duplicates in matrices for each dataset\n")
cat("#----------------------------------\n")
    
#create tmp dir to hold matrices
td <- tempdir()
print(td)   

g=4
#for( g in 1:4 ){

    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
        
    if( dir.exists(PATHS[2]) ){

        st1     = sprintf("%s/%s.csv.gz",PATHS[2],datasets[g,2])
        if( file.exists(st1) && file.info(st1)$size!=0 ){
            fin = data.table::fread(st1)
        }

        ids <- c()
        indx = which(colnames(fin)==IDS[ver,1])
        ids  = c(as.character(unlist(fin[,.SD,.SDcols=indx[1]])))
        dup  = duplicated(ids)
        map  = cbind(seq(1,length(ids),1),dup,ids)
                                    
            
        PATHS[5]=sprintf("%s/MATRIX",PATHS[2])
            
        cat("check for duplicates in ", datasets[g,1], " dataset...\n")
            
        ## can now read back matrix on disk from rds
        st1    = file.path(PATHS[5],sprintf("%s.rds",datasets[g,1]))
        mat    = readRDS(st1)
        NRold  = as.numeric(mat$nrow)
        NRnew  = as.numeric(table(map[,2])[[1]])

        if( NRnew != NRold ){

            cat("NR = ", NRold, " NR (unique) = ", NRnew, " correcting matrix size...\n")

            std = file.path(td,datasets[g,1])
            tmp = FBM(NRnew,NRnew,backingfile=std,FBMtype[FBMt])

            uni = as.numeric(map[map[,2]==FALSE,1])

            k=1
            for( i in 1:NRold ){
                if( map[i,2] == FALSE ){ 
                    tmp[,k] = mat[uni,i]
                    k=k+1
                }
            }

            ## save temp matrix
            tmp$save()
                
            ## clean our dataset space, if matrix already exist on disk.
            st1    = file.path(PATHS[5],sprintf("%s.rds",datasets[g,1]))
            if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }  
            st1    = file.path(PATHS[5],sprintf("%s.bk",datasets[g,1]))
            if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }  
                       

            ## copy temp matrix into our dataset space
            file.copy(from = file.path(td,sprintf("%s.rds",datasets[g,1])),
                      to   = file.path(PATHS[5],sprintf("%s.rds",datasets[g,1])))
            
            file.copy(from = file.path(td,sprintf("%s.bk",datasets[g,1])),
                      to   = file.path(PATHS[5],sprintf("%s.bk",datasets[g,1])))
   

            ## clean up
            unlink(sprintf("%s/%s",file.path(td),paste0(datasets[g,1] ,c(".bk",".rds")))) 
                
        }

        ## check that the path to the backing file points to correct location
        tmp = checkBK( tmp, file.path(PATHS[5],sprintf("%s.bk",datasets[g,1])) )

        ## save temp matrix
        tmp$save()
        
        ## clean up
        rm(mat)

        cat("...done.\n")
            
    }
#}#g
