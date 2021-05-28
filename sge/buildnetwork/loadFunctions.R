#----------------------------------
# load utility functions
#----------------------------------

#-- return average between two vectors 'x' and 'y'
##avg <- function(x,y) ((x+y)/2)
avg <- function(x,y){ ((x+y)/2) }

checkBK <- function(MAT, NEWbk){

    OLDbk = MAT$backingfile
    
    if( NEWbk != OLDbk ){ MAT$backingfile = NEWbk; }

    return(MAT)
    
}

unKnown <- function(x, DES="unknown"){
    x = ifelse(is.na(x),DES,x)
    x = ifelse(x=="",DES,x)
    return(x)
}

strMerge <- function(x){
    x = unique(x);
    x = ifelse(is.na(x),"",x);
    x = x[x!=""]

    if( length(x) == 0 ){ return(""); }
    if( length(x) == 1 ){ return(x);  }

    return(paste(x,collapse=";"));
    
}

findRACE <- function( FIN, find=c("hispanic"), targets=c("race_african_amer","race_white"), remove=c("race_more_than_one_calc")){
    
    ii = grep("race_",colnames(FIN));

    if( !is.null(find) ){
        ii = c(ii,grep(find,colnames(FIN)));
    }

    if( !is.null(remove) ){       
        rr = grep(remove,colnames(FIN));
        ii = ii[ifelse(is.na(match(ii,rr)),TRUE,FALSE)];
    }

    if( !is.null(targets) ){
        tt = match(targets,colnames(FIN));
        ii = ii[ifelse(is.na(match(ii,tt)),TRUE,FALSE)];
    }
        
    NR = dim(FIN)[1];
    NC = length(ii);
    
    xx = matrix(NA,nrow=NR,ncol=NC)
    colnames(xx) = colnames(FIN)[ii];

    xx2 = matrix(NA,nrow=NR,ncol=(length(targets)+1))
    colnames(xx2) = c(targets,"other")
    
    for( i in 1:NC ){
        xx[,i] = as.vector(unlist(FIN[,.SD,.SDcol=ii[i]]));
        #xx[,i] = ifelse(is.na(xx[,i]),"unknown",colnames(xx)[i])
    }

    x = rowSums(xx,na.rm=TRUE)
    x = ifelse(x==0,"",colnames(xx2)[(length(targets)+1)])
    
    for( i in 1:length(targets) ){
        xx2[,i] = as.vector(unlist(FIN[,.SD,.SDcol=tt[i]]));
        xx2[,i] = ifelse(is.na(xx2[,i]),"",colnames(xx2)[i])
    }

    xx2[,(length(targets)+1)] = x

    races = as.vector(apply(xx2,1,strMerge));
    races = ifelse(races=="","unknown",races);
    
    return(races);
    
}

#------------------------ NETWORK FUNCTIONS ------------------------

#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d = 1
    CC= as.numeric(length(V(dec[[1]])))
    for( i in 1:length(dec) ){
        if(as.numeric(length(V(dec[[i]]))) > CC){
            d=i
            CC=as.numeric(length(V(dec[[i]])))
        }
    }   

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

removeEdgeTerm <- function(GG,NAME){

    if( !is.null(igraph::get.edge.attribute(GG,NAME)) ){
        GG <- igraph::remove.edge.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.edge.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.edge.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

buildGraph <- function(VERTEX.ATTS, EDGE.ATTS){
## https://www.r-bloggers.com/2012/12/loading-huge-graphs-with-igraph-and-r/
    vertex.attrs = list(name=VERTEX.ATTS$name,
                        id=VERTEX.ATTS$id,
                        SPARKID=VERTEX.ATTS$SPARKID)

    edge.attrs = list(from=EDGE.ATTS$from,
                      to=EDGE.ATTS$to,
                      weight=EDGE.ATTS$weight,
                      zs=EDGE.ATTS$zs,
                      scaled=EDGES.ATTS$scales)

    edlist = rbind(match(EDGE.ATTS$from,vertex.attrs$name),
                   match(EDGE.ATTS$to,vertex.attrs$name))
    
    GG <- graph.empty(n=0,directed=F)
    GG <- add.vertices(GG,length(vertex.attrs$name),attr=vertex.attrs)
    GG <- add.edges(GG,edlist,attr=edge.attrs)

    rm(vertex.attrs,edge.attrs,edlist); gc();

    return(GG)

}

printGraphInfo <- function(GG,NAME){

    ##the number of ndoes, edges and density of graph
    N = as.numeric(length(V(GG)))
    E = as.numeric(length(E(GG)))

    ##print network's density
    rho = E/choose(N,2)
    cat("GRAPH = ", NAME, ": \n")
    cat("N = ",N, ", E = ", E, ", rho = ", rho,"\n")
    cat("--------\n")
    
}


setCnmax <- function( N, PER=10 ){

    #return( floor( (PER*N)/100) );
    return( floor(N/sqrt(N)) )
    
}


setCnmin <- function(  N ){ return(1) }

#---Get all edges internal to a community
intraEdges <- function(GG, ALG, CC, INTRA=NULL, INTER=NULL, useWE=NULL){

    intraINDX=NULL #edges in the community CC
    interINDX=NULL #edges going out from community CC   

    intra = NULL 
    inter = NULL 
    
    if( !is.null(igraph::get.vertex.attribute(GG,ALG)) ){

        coms <- get.vertex.attribute(GG,ALG)

        if( length(which(coms == CC)) != 0 ){
        
            ed_cc = E(GG)[inc(coms == CC)]

            all_edges_m <- get.edges(GG, ed_cc) #matrix representation           
            
            interINDX = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

            intraINDX = (ed_cc[(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])
            
        }
    }

    if( INTRA==TRUE && !is.null(intraINDX) ){
        intra_m = get.edges(GG,intraINDX)
        intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
        if( !is.null(useWE) ){
            intra = cbind(intra, get.edge.attribute(GG,useWE,intraINDX))
        }

        return(intra)
    }

    if( INTER==TRUE && !is.null(interINDX) ){
        inter_m = get.edges(GG,interINDX)        
        inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
        if( !is.null(useWE) ){
            inter = cbind(inter, get.edge.attribute(GG,useWE,interINDX))
        }

        return(inter)
    }

    return(NULL)
        
}


reCluster <- function( GG, ALGN, CnMAX, CnMIN, useWE=NULL ){

    if( !is.null(igraph::get.vertex.attribute(GG,ALGN)) ){

        #--- algorithm clustering 1
        ALG  <- get.vertex.attribute(GG,ALGN,V(GG))
        ALG1 <- cbind(V(GG)$name, ALG)

        Cn <- table(as.numeric(ALG1[,2]))    
        cc <- names(Cn)[Cn > CnMAX]


        RES <- list()
        k=1
        for( i in 1:length(cc) ){

            edCC = NULL
            edCC = intraEdges(GG, ALGN, cc[i], INTRA=TRUE, useWE=useWE )

            if( !is.null(edCC) ){

                ##ggLCC    <- graph_from_data_frame(d=edCC[,1:2], directed=F)

                # test using '_from_edgelist' instead of '_from_data_frame'
                ggLCC    <- graph_from_edgelist(el=edCC[,1:2], directed=F)

                #if( !is.null(useWE) ){
                #    ggLCC <- set.edge.attribute(ggLCC, useWE, E(ggLCC), as.numeric(edCC[,3]))
                #}
                
                #fc
                if( grepl("fc",ALGN) ){
                    if( !is.null(useWE) ){
                        res      <- igraph::cluster_fast_greedy(ggLCC, weights=as.numeric(edCC[,3]))         
                    } else {
                        res      <- igraph::cluster_fast_greedy(ggLCC, weights=NULL)
                    }
                    oo       <- cbind(res$names, res$membership)         
                }
                
                
                #louvain
                if( grepl("louvain",ALGN) ){
                    if( !is.null(useWE) ){
                        res      <- igraph::cluster_louvain(ggLCC, weights=as.numeric(edCC[,3]))
                    } else {                
                        res      <- igraph::cluster_louvain(ggLCC, weights=useWE)
                    }
                    oo       <- cbind(res$names, res$membership)         
                }

               
                RES[[k]]      <- oo
                names(RES)[k] <- cc[i]
                k=k+1
            }
            
        }#for


        if( length(RES) == 0 ){ return(NULL) }
        
        #--- algorithm clustering 2
        ALG2     <- cbind(ALG1, rep(-1, length(ALG1[,1])))
        indx     <- match(ALG2[,2],cc)
        indx     <- ifelse(is.na(indx),TRUE, FALSE)
        ALG2[,3] <- ifelse(indx, ALG2[,2], ALG2[,3])

        CCmax = max(as.numeric(ALG2[,3]))

        for( i in 1:length(RES) ){
 
            temp     <- RES[[i]]
            temp[,2] <- as.numeric(temp[,2]) + CCmax
        
            indx <- match(ALG2[,1],temp[,1])
            indx <- temp[indx,2]
        
            ALG2[,3] = ifelse(is.na(indx),ALG2[,3],indx)

            CCmax = max(as.numeric(ALG2[,3]))
        
        }

        #---reorder ALG2[,3]
        N = length(V(GG));
    
        temp    <- rep(-1, N)
        counter <- min(as.numeric(ALG2[,3]))
        Knew    <- 1;
        Kmax    <- max(as.numeric(ALG2[,3]))

        while( counter <= Kmax ){

            found=FALSE;
        
            for(v in 1:N ){
                if( as.numeric(ALG2[v,3]) == counter ){
                    temp[v] = Knew;
                    found=TRUE;
                }
            }

            if(found) Knew=Knew+1;
      
            counter=counter+1;
        }
        

        #---final 
        ALG3 <- cbind(ALG2, temp)
        return(ALG3)
    }

    return(NULL)
    
}

cluster2Level <- function( GG, ALG, CnMIN, CnMAX, nLEVELS=2, useWE=NULL ){
    
    for( levels in 1:nLEVELS ){

        if( levels == 1 ){
            ALGN <- ALG
        } else {
            ALGN <- sprintf("%s%d",ALG,levels)
        }

        levels <- levels + 1

        for( a in 1:length(ALGN) ){

            cat("reclustering ", ALGN[a], "...")

            oo = reCluster( GG, ALGN[a], CnMAX, CnMIN, useWE )

            if( !is.null(oo) ){
    
                #--- reclustering name
                ALGN2 <- sprintf("%s%d",ALG[a],levels)

                GG = removeVertexTerm(GG,ALGN2)
    
                if( is.null(igraph::get.vertex.attribute(GG,ALGN2)) ){
                    GG <- igraph::set.vertex.attribute(GG,ALGN2,V(GG),as.numeric(oo[,4]))
                }       
                
            }
        
            cat("...done.\n")
    
        }

    }#levels

    return(GG)
    
}

#---------------------------------------------------------------------------------

#------------------------ PLOT FUNCTIONS ------------------------

plotDist <- function( DF, XLAB="", YLAB="count", TIT="", DIR="./",PLOTn="XX"){

    DF = as.data.frame(DF)
    colnames(DF) = c("X","group")
    
    gplot <- ggplot(DF, aes(as.numeric(as.vector(X)), fill = group)) + geom_bar(show.legend=F)+
    labs(x=XLAB,y=YLAB,title=TIT)+        
    scale_fill_manual(values = c("LOWER" = "#FA42DD",
                                 "NON"   = "#BEBBC7",
                                 "UPPER" = "#5039FA"))+
        theme(
                title=element_text(face="bold",size=rel(1.5)),
                axis.title.x=element_text(face="bold",size=rel(2.0)),
                axis.title.y=element_text(face="bold",size=rel(2.0)),
                legend.title=element_text(face="bold",size=rel(2.0)),
                legend.text=element_text(face="bold",size=rel(2.0)),
                legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_rect(linetype="solid",fill=NA))


    #---WIDTH and HEIGHT for plots
    WIDTH=480
    HEIGHT=480

    png(sprintf("%s/%s_dist.png",DIR,PLOTn),width=WIDTH,height=HEIGHT,units="px")
    print(gplot)
    dev.off()

    return(gplot)
}

plotNN <- function( DF, XLAB="", YLAB="count", DIR="./", PLOTn="XX"){

     xpoint = max(as.numeric(as.vector(DF[DF[,2]=="NN",1])))
        
    gplot <- ggplot(DF, aes(as.numeric(as.vector(DF[,1])),fill=DF[,2])) +
        geom_histogram(bins=20, colour="black",show.legend=F)+
    labs(x=XLAB,y=YLAB,title="")+        
        scale_fill_manual(values = c("NN"    = "#FA42DD",
                                     "NON"   = "#BEBBC7"))+
        geom_vline(xintercept=xpoint, linetype="dashed", color = "#FA42DD", size=rel(1.5))+
        theme(title=element_text(face="bold",size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(2.0)),
              legend.text=element_text(face="bold",size=rel(2.0)),
              legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_rect(linetype="solid",fill=NA))


    #---WIDTH and HEIGHT for plots
    WIDTH=480
    HEIGHT=480

    png(sprintf("%s/%s_degree.png",DIR,PLOTn),width=WIDTH,height=HEIGHT,units="px")
    print(gplot)
    dev.off()

}




#---------------------------------------------------------------------------------

