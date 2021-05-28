
## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

run <- vector(length=2)
run[1] = 1 # build adjacency matrices for each dataset
run[2] = 0 # build network from joint dataset
run[3] = 0 # build fused network 

if( run[1] ){

    #----------------------------------
    # build adjacency matrices for each dataset
    #----------------------------------

    cat("#----------------------------------\n")
    cat("# build adjacency matrices for each dataset\n")
    cat("#----------------------------------\n")    

    for( g in 1:4 ){

        PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[g,1])
        
        if( dir.exists(PATHS[2]) ){

            PATHS[5]=sprintf("%s/MATRIX",PATHS[2])
            
            cat("build adjacency matrix for ", datasets[g,1], " dataset...\n")

            ## can now read back matrix on disk from rds
            mat = readRDS(file.path(PATHS[5],sprintf("%s.rds",datasets[g,1])))       
            NR  = as.numeric(mat$nrow)
            #---done

            ## temp code
            ## check that the path to the backing file points to correct location
            mat = checkBK( mat, file.path(PATHS[5],sprintf("%s.bk",datasets[g,1])) )
            
            ## center and log each patient's similarity scores, i.e.
            ## to obtain Gaussian distribution 
            ii=seq_along(mat[,1])
            Mn=rep(NA,length(ii))
            Sd=rep(NA,length(ii))

            for( i in 1:length(ii) ){
                x     = mat[,i]
                x     = log(x)
                x     = x[!is.infinite(x)]
                Mn[i] = mean(x,na.rm=T)
                Sd[i] = sd(x,na.rm=T)
            }
            #---done

            #create our adjacency matrix 
            td <- PATHS[5]
            print(td)
            st1 = sprintf("%s/%s_ANN.bk",td,datasets[g,1])
            if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }
            st1 = sprintf("%s/%s_ANN.rds",td,datasets[g,1])
            if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }

            ann = FBM(NR,NR,backingfile=file.path(td,sprintf("%s_ANN",datasets[g,1])), FBMtype[FBMt])
            #---done


            # select the nearest neighbours for each patient, by assuming each patients
            # similarity score is gaussian distributed,
            # and then apply a cut-off (p-value) to the tail of the distribution. 
            for( i in 1:length(ii) ){
                x    = log(mat[,i])
                zs   = (x-Mn[i])/Sd[i]
                pval = 2*pnorm(-abs(zs)) 
                ann[,i] = 0
                ann[i,] = 0
                ann[ii[pval <= alpha[thres] & pval>0.0 & zs <= 0],i] = 1
                ann[i,] = ann[,i]
            }
            #---done


            ## create rds files
            ann$save()        

        }#if
            
    }#g

    cat("#----------------------------------\n")
    cat("# done                             \n")
    cat("#----------------------------------\n")
    
}

if( run[2] ){

    #----------------------------------
    # build network for dataset
    #----------------------------------

    cat("#----------------------------------\n")
    cat("# build network for dataset\n")
    cat("#----------------------------------\n") 

    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[choice,1])
    
    if( dir.exists(PATHS[2]) ){

        PATHS[5]=sprintf("%s/MATRIX",PATHS[2])
        PATHS[6]=sprintf("%s/%s",PATHS[3],"GRAPHS")
        
        cat("build network for ", datasets[choice,1], " dataset...\n")   
    
        ## can now read back matrix on disk from rds
        mat = readRDS(file.path(PATHS[5],sprintf("%s.rds",datasets[choice,1])))
        NR  = as.numeric(mat$nrow)
        #---done

    ## center and log each patient's similarity scores, i.e. to obtain Gaussian distribution 
        ii=seq_along(mat[,1])
        Mn=rep(NA,length(ii))
        Sd=rep(NA,length(ii))
        mx=rep(NA,length(ii))
        mn=rep(NA,length(ii))

        for( i in 1:length(ii) ){
            x     = mat[,i]
            mx[i] = max(x,na.rm=T)
            mn[i] = min(x,na.rm=T)

            x     = log(x)
            x     = x[!is.infinite(x)]
            Mn[i] = mean(x,na.rm=T)
            Sd[i] = sd(x,na.rm=T)
        }
        #---done

        #create our adjacency matrix 
        td <- tempdir()
        print(td)
        adj = FBM(NR,NR,backingfile=file.path(td,"adj"),FBMtype[FBMt])
        #---done

    # select the nearest neighbours for each patient, by assuming each patients similarity score is gaussian distributed,
    # and then apply a cut-off (p-value) to the tail of the distribution. 
        for( i in 1:length(ii) ){
            x    = log(mat[,i])
            zs   = (x-Mn[i])/Sd[i]
            pval = 2*pnorm(-abs(zs))
            adj[,i] = 0
            adj[i,] = 0
            adj[ii[pval <= alpha[thres] & pval>0.0 & zs <= 0],i] = 1
            adj[i,] = adj[,i]
        }
        #---done



        ##create edge-list from adj (for igraph)
        df = data.frame(ii=as.numeric(),jj=as.numeric())
        for( i in 1:length(ii) ){
            dfj  = ii[adj[,i]==1]
            dfi  = rep(i,length=length(dfj))
            temp = data.frame(ii=dfi,jj=dfj)
            df   = rbind(df,temp)
        }
        #---done

        ## clean up
        unlink(sprintf("%s/%s",file.path(td),paste0("adj", c(".bk", ".rds"))))
        #---done

        ##building graph (using igraph)
        ###gg  = graph_from_data_frame(df[,1:2],directed=F)
        gg  = graph_from_edgelist(df[,1:2],directed=F)
        printGraphInfo(gg,"raw")
        #---done

        ## read in the dataset
        fin           <- list()
        fin[[1]]      <- data.table::fread(sprintf("%s/%s.csv.gz",PATHS[2],datasets[choice,2]))
        names(fin)[1] <- datasets[choice,1]
        #---done

        # SPARK patient id
        gg  = removeVertexTerm(gg,"SPARKID")
        id  = as.character(unlist(fin[[1]][,1]))
        ii  = seq(1,length(id),1)
        val = id[match(V(gg)$name,ii)]
        gg  = igraph::set.vertex.attribute(gg,"SPARKID",V(gg),val)
        #---done

        # SPARK family id
        gg  = removeVertexTerm(gg,"SPARKFAMILYID")
        fid = as.character(unlist(fin[[1]][,2]))
        ii  = seq(1,length(fid),1)
        val = fid[match(V(gg)$name,ii)]
        gg  = igraph::set.vertex.attribute(gg,"SPARKFAMILYID",V(gg),fid)
        #---done

        # simplify network
        gg2 = simplify(gg,remove.multiple=T,remove.loops=T)

        printGraphInfo(gg2,"simplified")

        ##write network to file (in gml format)
        write.graph(gg2,sprintf("%s/%s_%s.gml",PATHS[6],datasets[choice,1],aname[thres]),format="gml")
        #---done
    
        # find network's largest connected component
        gg2 = findLCC(gg2)

        printGraphInfo(gg2,"LCC")
        #---done
    
        # clustering, using louvain and fast-greedy algorithms
        lv = cluster_louvain(gg2)
        fc = cluster_fast_greedy(gg2)

        gg2=removeVertexTerm(gg2,"louvain")
        gg2=set.vertex.attribute(gg2,"louvain",V(gg2),lv$membership)
        gg2=removeVertexTerm(gg2,"fc")
        gg2=set.vertex.attribute(gg2,"fc",V(gg2),fc$membership)
        #---done

        ##write networks LCC to file (gml format)
        write.graph(gg2,sprintf("%s/%s_%s_LCC.gml",PATHS[6],datasets[choice,1],aname[thres]),format="gml")
        #---done

        cat("#----------------------------------\n")
        cat("# done                             \n")
        cat("#----------------------------------\n")

    }#if
        
}


if( run[3] ){


    #----------------------------------
    # build network for fused datasets
    #----------------------------------

    #----------------------------------
    # compile c++ code (fuse.cpp)
    #----------------------------------

    PATHS[2]=sprintf("%s/%s",PATHS[1],datasets[choice,1])
    
    if( dir.exists(PATHS[2]) ){
                
        #Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
        #Sys.setenv("PKG_LIBS" = "-fopenmp")
        #Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -O3 -DLEVEL1_DCACHE_LINESIZE=`getconf LEVEL1_DCACHE_LINESIZE` -I\"/usr/include/gperftools/\" -ltcmalloc_minimal")
        #Sys.setenv("PKG_LIBS" = "-fopenmp -L\"/usr/lib/x86_64-linux-gnu/\" -ltcmalloc_minimal")
        Sys.setenv("PKG_CXXFLAGS" = cxxflag.opts)
        Sys.setenv("PKG_LIBS"     = lib.opts)
        sourceCpp(sprintf("%s/%s",PATHS[4],files[6]))
    
        cat("#----------------------------------\n")
        cat("# build fused network\n")
        cat("#----------------------------------\n")

        PATHS[5]=sprintf("%s/MATRIX",PATHS[3])

        fuseM1 = readRDS(file.path(PATHS[5],sprintf("%s.rds",fuseName[1])))
        fRes1  = read.delim(file.path(PATHS[5],sprintf("%s.csv",fuseName[1])),sep=" ",header=T)

        cat("Get edges from fuse matrix...\n")
        
        # build network for fuseM1
        NR  = as.numeric(fuseM1$nrow)
        res = getEdgelist( fuseM1, nCORES, c(1) );
        df  = res$RESULT

        cat("done.\n")
        

        ## get all node ids in df, order from lowest to highest
        node_ids_df1 = unique(df[,1])
        node_ids_df2 = unique(df[,2])
        node_ids_df  = unique(c(node_ids_df1,node_ids_df2))
        node_ids_df  = node_ids_df[order(node_ids_df,decreasing=F)]
        
        ## read in spark IDs
        ids = read.delim(sprintf("%s/%s.csv",PATHS[3],files[8]),sep="\t", header=F)[[1]]
        ids = as.vector(ids)

        ## mapping for SPARK patient ids
        spark_ids = cbind(ids,seq_along(1:length(ids)))

        ## match SPARK patients ids, with node ids in df
        indx = match(node_ids_df,as.numeric(spark_ids[,2]))

        ## build node attribute data frame for df
        nodes_df = cbind(node_ids_df,seq_along(1:length(node_ids_df)),spark_ids[indx,1])
        nodes_df = as.data.frame(nodes_df)
        colnames(nodes_df) <- c("name","id","SPARKID")
        #----

        ##garbage collection nodes
        rm(node_ids_df1,node_ids_df2,node_ids_df,ids,spark_ids,indx); gc();
        
        ## build edge attribute data frame for df, and add edge weighting schemes (i.e. znorm and scaled)
        x  <- df[,3]
        zs <- (log(x) - mean(log(x)))/sd(log(x)) #z-normed             
        sc <- (zs - min(zs))/(max(zs)-min(zs))   #scaled

        df       = cbind(df,zs,sc)
        edges_df = as.data.frame(df)
        colnames(edges_df) <- c("from", "to", "weight","zs","scaled")
        ##------

        ##garbage collection edges
        rm(x,zs,sc,df); gc();
        
        # convergence plots
        PATHS[5]=sprintf("%s/PLOTS",PATHS[3])

        N=3:length(fRes1[,1])    
        png(filename=file.path(PATHS[5],sprintf("%s.png",fuseName[1])), width=WIDTH, height=HEIGHT, units="px")
        plot(xy.coords(fRes1[N,1],fRes1[N,4]),xlab="steps",ylab="Pc difference",main=sprintf("ANN (p <= %s cut-off applied to datasets)",aname[thres]),col="blue")
        dev.off()
       

        cat("build network in igraph...")
        
        # building graph (using igraph)
        ##gg  = graph_from_data_frame(df[,1:2],directed=F) #(note this is slow)
        ##gg  = graph_from_edgelist(as.matrix(df[,1:2,drop=FALSE],directed=F))
        ##gg = buildGraph(VERTEX.ATTS=nodes_df,EDGE.ATTS=edges_df)
        gg = graph_from_data_frame(d=edges_df, vertices=nodes_df, directed=FALSE)

        rm(nodes_df,edges_df); gc();
        
        cat("done.\n")
        
        printGraphInfo(gg,"raw")

        # add edge weights to graph
        #x <- df[,3]
        #removeEdgeTerm(gg,"weight")
        #gg <- set.edge.attribute(gg,"weight",E(gg),x)

        #zs <- (log(x) - mean(log(x)))/sd(log(x))
        #removeEdgeTerm(gg,"znorm")
        #gg <- set.edge.attribute(gg,"znorm",E(gg),zs)
    
        #sc <- (zs - min(zs))/(max(zs)-min(zs))
        #removeEdgeTerm(gg,"scaled")
        #gg <- set.edge.attribute(gg,"scaled",E(gg),sc)
    
        # read in spark IDs
        #ids = read.delim(sprintf("%s/%s.csv",PATHS[3],files[8]),sep="\t", header=F)[[1]]
        #ids = as.vector(ids)

        # SPARK patient id
        #ii  = seq_along(1:length(ids))
        #gg  = removeVertexTerm(gg,"SPARKID")    
        #val = as.character(ids[match(V(gg)$name,ii)])
        #gg  = igraph::set.vertex.attribute(gg,"SPARKID",V(gg),val)

        cat("cluster igraph network...\n")
        
        # clustering, using fast-greedy algorithms
        fc   = cluster_fast_greedy(gg, weights=NULL)
        gg=set.vertex.attribute(gg,"fc",V(gg),fc$membership)
        cat(" no weight.\n")

        fcWe = cluster_fast_greedy(gg, weights=E(gg)$weight)
        gg=set.vertex.attribute(gg,"fcWe",V(gg),fcWe$membership)
        cat(" with weight.\n")
        
        fcSc = cluster_fast_greedy(gg, weights=E(gg)$scaled)
        gg=set.vertex.attribute(gg,"fcSc",V(gg),fcSc$membership)
        cat(" with scaled weight.\n")
        
        #gg=removeVertexTerm(gg,"fc")
        #gg=removeVertexTerm(gg,"fcWe")
        #gg=removeVertexTerm(gg,"fcSc")


        cat("...done.\n")
    
        ## write network to file (in gml format)
        PATHS[6]=sprintf("%s/%s",PATHS[3],"GRAPHS")    
        write.graph(gg,sprintf("%s/%s.gml",PATHS[6],ggName[1]),format="gml")
     
    
        cat("#----------------------------------\n")
        cat("# done                             \n")
        cat("#----------------------------------\n")

    }#if
}
