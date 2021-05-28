#----------------------------------
# load variable selection functions
#----------------------------------

#--- set random number seed
setSeed <- function(SEED){
    if( is.null(SEED) ){
        SEED = as.integer(Sys.time())
        set.seed(SEED)
    } else {
        set.seed(SEED)
    }

    return(SEED)
}

Min <- function(x){ min(x,na.rm=T)  }
Max <- function(x){ max(x,na.rm=T)  }

Extract <- function(x, ii){ x[ii] }

dealWithNA <- function(MAT,NAval=-1){
    MAT[which(is.na(MAT),arr.ind=TRUE)]=NAval
    return(MAT)
}

rough.fixNA <- function(MAT){
## fill NA with it's column median
## Ref: https://cran.r-project.org/web/packages/RRF/RRF.pdf

    indx = which(is.na(MAT),arr.ind=TRUE)
    Cols = unique(indx[,2])
    for( i in 1:length(Cols) ){
        val = MAT[,Cols[i]]
        med = median(val,na.rm=T)
        MAT[,Cols[i]] = ifelse(is.na(val),med,val)
    }

    return(MAT)
    
}

fscore <- function(X){

    classes <- levels(as.factor(X[,1]))
    pred    <- levels(as.factor(X[,2]))

    Nclasses <- length(classes)
    Npred    <- length(pred)

    ## binary confusion matrix
    ##      asd   con
    ## asd  TP    FN
    ## con  FP    TN
    ##------------------------ 
    conmat  <- matrix(NA,nrow=Nclasses,ncol=Npred)
    rownames(conmat) <- classes
    colnames(conmat) <- pred

    for( i in 1:Nclasses ){
        for( j in 1:Npred ){
            conmat[i,j] = sum(X[,1] == classes[i] & X[,2] == pred[j])
        }
    }

    ## precision = TP / (TP+FP)
    prec        = rep(NA,(Nclasses+1))
    names(prec) = c(classes,"total")

    ## recall = TP / (TP+FN)
    recall = rep(NA,Nclasses+1)
    names(recall) = c(classes,"total")
    
    TP  = diag(conmat)   
    
    prec[1] = TP[1]/sum(conmat[,1]);
    prec[2] = TP[2]/sum(conmat[2,]);
    prec[3] = sum(TP)/sum(conmat);

    recall[1] = TP[1]/sum(conmat[1,]);
    recall[2] = TP[2]/sum(conmat[,2]);
    recall[3] = sum(TP)/sum(conmat);

    Fscore = (2*prec*recall)/(prec+recall);
    names(Fscore) = c(classes,"total")
    
    return(list(conmat=conmat,prec=prec,recall=recall,Fscore=Fscore))
}

logBase <- function(x, BASE=0){

    x = as.numeric(x)

    if( BASE == 2 ){
        return(log2(x))
    } else if ( BASE == 10 ){
        return(log10(x))
    } 

    return(log(x))

}

lgrad <- function(n,k,p){
    ## Ref: https://mc-stan.org/docs/2_21/functions-reference/binomial-distribution.html
    return( k/(p+SMALL) - (n-k)/(1-p+SMALL)  );
}

lbin <- function(n,k){
    ## Ref: https://mc-stan.org/docs/2_21/functions-reference/binomial-distribution.html
    return(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1))
}

lbinomial <- function(n,k,p){
    ## Ref: https://mc-stan.org/docs/2_21/functions-reference/binomial-distribution.html
    return( lbin(n,k) + k * log(p+SMALL) + (n-k) * log(1-p+SMALL) );
}

milog <- function(pxyz, pxz, pyz, pz, BASE=0){

    #SMALL <- 1.0e-100
    
    if( pxyz == 0 ){ return(0) }

    if( pxz == 0 ){ pxz = pxz + SMALL }
    if( pyz == 0 ){ pyz = pyz + SMALL }
    if( pz  == 0 ){ pz  = pz  + SMALL }
    
    return(pxyz*logBase( (pxyz*pz)/(pxz*pyz), BASE ))
    
}

entropy <- function(x, BASE=0){
    x   <- as.numeric(x)    
    if( is.na(x) ){ return(NA) }
    else if( x == 0 ){ return(0) }
    else { return(-x*logBase(x,BASE)) }    
}

MutualInformation <- function(X,Y,BASE=0){

    RES <- list()

    Hx    = 0
    Hy    = 0
    Hxy   = 0
    MIxy  = 0
    NMIxy = 0
    
    X = as.numeric(X)
    Y = as.numeric(Y)
    
    indx = !is.na(X) & !is.na(Y)

    x = X[indx]
    y = Y[indx]

    if( length(x) > 0 && length(y) > 0 ){
    
        ct2 = table(x,y)        
        N   = sum(ct2)
        ix  = nrow(ct2)
        jy  = ncol(ct2)

        # H(X)
        for( i in 1:ix ){
            nx = sum(ct2[i,])
            px = nx/N
            Hx = Hx + entropy(px,BASE)
        }       
        
        # H(Y) and H(X|Y)
        for( j in 1:jy ){
            ny = sum(ct2[,j])
            py = ny/N
            Hy = Hy + entropy(py,BASE)            
            sum = 0
            for( i in 1:ix ){
                nxy = ct2[i,j]
                pxy = nxy/ny
                sum = sum + entropy(pxy,BASE)
            }            
            Hxy = Hxy + py*sum
        }

        MIxy  = Hx - Hxy
        NMIxy = 2*MIxy/(Hy+Hx)
        
    }
        
    RES$Ha  = Hx
    RES$Hb  = Hy
    RES$Hab = Hxy
    RES$MI  = MIxy
    RES$NMI = NMIxy
    
    return(RES)
    
}

conditionalMI <- function(X,Y,Z,BASE=0){

    RES <- list()

    MIxyCz = 0
    MIxzCy = 0
    MIxyz  = 0
    IIxyz  = 0
    SRxyz  = 0
    
    X = as.numeric(X)
    Y = as.numeric(Y)
    Z = as.numeric(Z)    
    
    indx = !is.na(X) & !is.na(Y) & !is.na(Z)

    x = X[indx]
    y = Y[indx]
    z = Z[indx]

    if( length(x) > 0 && length(y) > 0 && length(z) > 0 ){

        MIxy = MutualInformation(X,Y,BASE)
        MIxz = MutualInformation(X,Z,BASE)
        MIyz = MutualInformation(Y,Z,BASE)
                
        ctxyz = table(x,y,z)        
        N     = sum(ctxyz)
        II    = dim(ctxyz)

       
        #-----------------------------------------
        # Interaction Information:  II(X;Y;Z) = I(X;Y) - I(X;Y|Z)
        # II(X;Y;Z) = I(X;Y) - I(X;Y|Z)
        #-----------------------------------------
        # Conditional MI:            I(X;Y|Z) (is symmetric)
        # I(X;Y|Z) = H(X;Z) + H(Y;Z) - H(X,Y,Z) - H(Z)
        #          = H(X|Z) + H(Y|Z) - H(X,Y,Z) + H(Z)
        #          = Sum_{x,y,z} p(x,y,z) * log( (p(x,y,z)*p(z))/(p(x,z)*p(y,z)) )
        #-----------------------------------------
        # Joint Mutual Information:  I(X,Y;Z) 
        # I(X,Y;Z) = I(X;Z) + I(Y;Z) + I(X;Y|Z) - I(X;Y)
        #-----------------------------------------
        # Refs: https://math.stackexchange.com/questions/168316/calculating-conditional-entropy-given-two-random-variables
        #     : https://en.wikipedia.org/wiki/Mutual_information
        #     : https://en.wikipedia.org/wiki/Interaction_information
        #     : https://www.mdpi.com/1099-4300/21/9/869/htm
        #     : https://arxiv.org/pdf/2001.06089.pdf
        #-----------------------------------------

        MIxyCz = 0
        Hxyz   = 0
        Hyz    = 0
        Hxz    = 0
        Hz     = 0        
        
        for( zi in 1:II[3] ){
            nz = sum(ctxyz[,,zi])
            pz = nz/N
            Hz = Hz + entropy(pz,BASE)
            for( yi in 1:II[2] ){
                nyz = sum(ctxyz[,yi,zi])
                pyz = nyz/N               
                Hyz = Hyz + entropy(pyz,BASE)                
                for( xi in 1:II[1] ){                                                      
                    nxyz = ctxyz[xi,yi,zi]
                    pxyz = nxyz/N                    
                    Hxyz = Hxyz + entropy(pxyz,BASE)                   
                }
            }
            for( xi in 1:II[1] ){
                nxz  = sum(ctxyz[xi,,zi])
                pxz  = nxz/N
                Hxz = Hxz + entropy(pxz,BASE)  
            }
        }        

        #Hy=0
        #Hxy=0
        #for( yi in 1:II[2] ){
        #    ny = sum(ctxyz[,yi,])
        #    py = ny/N
        #    Hy = Hy + entropy(py,BASE)
        #    for( xi in 1:II[1] ){
        #        nxy = sum(ctxyz[xi,yi,])
        #        pxy = nxy/N
        #        Hxy = Hxy + entropy(pxy,BASE)
        #    }
        #}
        
        # Conditional MI 
        #MIxyCz = MIxz$Hab + MIyz$Hab - Hxyz + MIxz$Hb#Hz
        MIxyCz = Hxz + Hyz - Hxyz - Hz

        # Correct for any rounding error
        if( MIxyCz < 0 ){ MIxyCz = 0 }
        
        # Joint MI
        MIxyz  = MIxz$MI + MIyz$MI + MIxyCz - MIxy$MI
        ##MIxyz2 = Hxy - Hxyz + Hz
        
        # Interaction Information
        IIxyz = MIxy$MI - MIxyCz                

        # SR = I(X,Y;Z)/H(X,Y,Z)
        SRxyz = MIxyz/Hxyz
        
        # Normalised Conditonal Mutual Information
        # I(X,Y|Z)/H(Y|Z)
        NMI = MIxyCz/Hyz
    }

    RES$MIabCc  <- MIxyCz   
    RES$MIacb   <- MIxyz
    #RES$MIacb2  <- MIxyz2
    RES$IIabc   <- IIxyz
    RES$Habc    <- Hxyz
    RES$SRabc   <- SRxyz
    RES$NMI     <- NMI
    RES$MIxy    <- MIxy
    RES$MIxz    <- MIxz
    RES$MIyz    <- MIyz
    #RES$HxCz    <- MIxz$Hab
    #RES$HyCz    <- MIyz$Hab
    #RES$Hz1     <- MIxz$Hb
    #RES$Hz2     <- MIyz$Hb
    #RES$Hx      <- MIxz$Ha
    #RES$Hy      <- MIyz$Ha
    #RES$Hz      <- Hz
    #RES$Hxz     <- Hxz
    #RES$Hyz     <- Hyz
    
    return(RES)
    
}


## feature selection
Jcmi <- function( S, IXkY, IXjXk, IXjXkCY, beta=1, gamma=1 ){

    Nf   <- length(S[,1])
    Sin  <- seq_along(1:Nf)[as.logical(S[,2])]
    Sout <- seq_along(1:Nf)[as.logical(S[,3])]    

    Jcmi <- rep(NA,Nf)

    if( length(Sin) > 0 && length(Sout) > 0 ){
    
    for( k in 1:length(Sout) ){
        val1  = as.numeric(IXkY[Sout[k]])
        val1  = ifelse( is.na(val1), 0, val1)
        term1 = val1#as.numeric(IXkY[Sout[k]])
        term2 = 0
        term3 = 0
        for( j in 1:length(Sin) ){
            val2 = IXjXk  [Sin[j],Sout[k]]
            val2 = ifelse( is.na(val2), 0, val2)

            val3 = IXjXkCY[Sin[j],Sout[k]]
            val3 = ifelse( is.na(val3), 0, val3)
            
            term2 = term2 + val2#IXjXk  [Sin[j],Sout[k]]
            term3 = term3 + val3#IXjXkCY[Sin[j],Sout[k]]
        }
        Jcmi[Sout[k]] = term1 - (beta * term2) + (gamma * term3)
    }

    }
        
    return(Jcmi)
}

## feature selection
Fcbf <- function( S, IXkY, IXjXk, beta=1 ){

    Nf   <- length(S[,1])
    Sin  <- seq_along(1:Nf)[as.logical(S[,2])]
    Sout <- seq_along(1:Nf)[as.logical(S[,3])]    

    fcbf <- rep(NA,Nf)

    if( length(Sin) > 0 && length(Sout) > 0 ){
    
    for( k in 1:length(Sout) ){
        val1  = as.numeric(IXkY[Sout[k]])
        val1  = ifelse( is.na(val1), 0, val1)
        term1 = val1#as.numeric(IXkY[Sout[k]])
        term2 = 0
        
        for( j in 1:length(Sin) ){
            val2 = IXjXk  [Sin[j],Sout[k]]
            val2 = ifelse( is.na(val2), 0, val2)                       
            term2 = term2 + val2
           
        }
        fcbf[Sout[k]] = term1 - (beta * term2)
    }

    }
        
    return(fcbf)
}



jointEntropy <- function(x,y){

    x    <- as.numeric(x)
    y    <- as.numeric(y)
    indx <- !is.na(x) & !is.na(y)
    x    <- x[indx]
    y    <- y[indx]
    N    <- length(x)
    sum  <- 0
    for( i in 1:N ){
        for( j in 1:N ){
            if( i <= j ){
                x.y <- x[i]*y[j]
                if( !is.na(x.y) ){
                    if( i == j ){
                        sum <- sum + entropy(x.y)
                    } else {
                        sum <- sum + 2*entropy(x.y)
                    }
                }
            }
        }
    }
    
    return(sum)
}


printPer6Months <- function(VATTS, DES, TARGET,PRINT=TRUE){

    N = length(VATTS[,1])

    CN = colnames(VATTS)

    Indx = match(TARGET,CN)

    val = as.vector(VATTS[,Indx[1]])
    val = as.numeric(ifelse(val=="unknown",-1,val))

    perMonths = c("unknown",6,12,18,24,30,36,42,48,52)
    months = matrix(NA,ncol=length(perMonths),nrow=1)
    colnames(months) = perMonths
    
    x = rep("",N)

    for( i in 1:length(perMonths) ){

        Mrange = colnames(months)[i]
        tally  = 0 
        
        switch(Mrange,
               "unknown" = { indx    = val < 0;    },
               "6" ={ indx = val > 0  & val <= 6;  },               
               "12"={ indx = val > 6  & val <= 12; },
               "18"={ indx = val > 12 & val <= 18; },
               "24"={ indx = val > 18 & val <= 24; },
               "30"={ indx = val > 24 & val <= 30; },
               "36"={ indx = val > 30 & val <= 36; },
               "42"={ indx = val > 36 & val <= 42; },
               "48"={ indx = val > 42 & val <= 48; },
               "52"={ indx = val > 52;             }
               )

        x[indx]     = Mrange;        tally       = sum(indx);
        months[1,i] = tally;
        
    }
    
    if( PRINT ){
    cat("---------------\n")
    cat("", DES, " variable: ", colnames(VATTS)[Indx[1]], "\n" );
    print(months)
    months = (months/N)*100
    print(months)    
    #print((table(months)/N)*100)
    cat("---------------\n")
    }
    
    return(x);
    
}

printVattVariables <- function(VATTS, DES, TARGETS){

    N = length(VATTS[,1])

    CN = colnames(VATTS)

    indx = match(TARGETS,CN)

    if( length(indx) != 0 ){
        for( i in 1:length(indx) ){
            cat("", DES, " variable: ", colnames(VATTS)[indx[i]], "\n" );
            print(table(VATTS[,indx[i]]))
            print((table(VATTS[,indx[i]])/N)*100)
            cat("----- DONE -----\n")
        }
    }
    
}

buildVattsTable <- function(GG, VATTS){

    NR = length(V(GG))
    NC = length(VATTS)
    oo=matrix(NA,NR,NC)
    colnames(oo) = VATTS

    for( i in 1:NC ){
        x = igraph::get.vertex.attribute(GG,VATTS[i])
        if( !is.null(x) ){
            oo[,i] = x
        }
    }

    return(oo)
    
}

addPrior <- function(x, alpha){
    x[x==0] = alpha
    x = x/sum(x)
}

shuffle <- function(MAT,DIAG=FALSE,P=1){
    ## randomly shuffle P% elements in upper and lower diagonal MAT 
    X = NULL
    
    if( !is.null(MAT) ){
        NR = dim(MAT)[1]
        NC = dim(MAT)[2]

        if( P==0 ){ return(MAT); }
        
        if( NR == NC && NR > 1 ){
        
            if( P == 1 ){ #shuffle all
                indx      = upper.tri(MAT,diag=DIAG)
                MAT[indx] = sample(MAT[indx])
                indx      = lower.tri(MAT,diag=DIAG)
                MAT[indx] = sample(MAT[indx])
                X         = MAT
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            } else {#shuffle fraction (P) of elements
                indx              = upper.tri(MAT,diag=DIAG)
                ele               = cbind(which(indx,arr.ind=T),MAT[indx])
                Nele              = dim(ele)[1]
                Subele            = ele[sample(Nele,floor(Nele*P),replace=F),]
                Subele[,3]        = Subele[sample(dim(Subele)[1]),3]
                MAT[Subele[,1:2]] = Subele[,3]
                indx              = lower.tri(MAT,diag=DIAG)
                ele               = cbind(which(indx,arr.ind=T),MAT[indx])
                Nele              = dim(ele)[1]
                Subele            = ele[sample(Nele,floor(Nele*P),replace=F),]
                Subele[,3]        = Subele[sample(dim(Subele)[1]),3]
                MAT[Subele[,1:2]] = Subele[,3]
                X                 = MAT
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            }
        }
    }

    return(X)
    
}

shuffle2 <- function(MAT,P=1){
    ## randomly shuffle P% elements in MAT 
    X = NULL

    if( !is.null(MAT) ){
        NR = dim(MAT)[1]
        NC = dim(MAT)[2]

        if( P == 0 ){ return(MAT); }
        
        if( P == 1 ){
        
            if( NR == 1 ){          
                X = MAT[1,sample.int(NC)]
                names(X) = colnames(MAT)          
            } else {
                X  = MAT[sample.int(NR),sample.int(NC)]
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            }

        } else {

            if( NR == 1 ){
                X          = matrix(NA,nrow=NR,ncol=NC)
                Indx       = sample.int(floor(P*NC))
                RIndx      = sample(Indx,replace=F)
                X[1,RIndx] = MAT[1,Indx] 
                X[1,]      = ifelse(is.na(X[1,]),MAT[1,],X[1,])
                #names(X)   = colnames(MAT)
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            } else {
                X       = matrix(NA,nrow=NR,ncol=NC)
                Indx    = seq(1,(NR*NC),1)
                Indx    = sample(Indx, floor(P*NR*NC),replace=F)
                Rndx    = sample(Indx)
                X[Rndx] = MAT[Indx]
                indx    = which(is.na(X),arr.ind=T)
                X[indx] = MAT[indx]
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            }
        }#ifelse           
    }##null
    
    return(X)

}

shuffle3 <- function(MAT,DIAG=FALSE,P=1){
    ## randomly select P% elements in upper and lower diagonal MAT and then set values randonly between [0,1]
    X = NULL
    
    if( !is.null(MAT) ){
        NR = dim(MAT)[1]
        NC = dim(MAT)[2]

        if( P == 0 ){ return(MAT); }
        
        if( NR == NC && NR > 1 ){

            if( P == 1 ){
                indx      = upper.tri(MAT,diag=DIAG)
                MAT[indx] = runif(length(indx),0,1)
                indx      = lower.tri(MAT,diag=DIAG)
                MAT[indx] = runif(length(indx),0,1)
                X         = MAT
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            } else {
                indx              = upper.tri(MAT,diag=DIAG)
                ele               = cbind(which(indx,arr.ind=T),MAT[indx])
                Nele              = dim(ele)[1]
                Subele            = ele[sample(Nele,floor(Nele*P),replace=F),]
                Subele[,3]        = runif(dim(Subele)[1],0,1)
                MAT[Subele[,1:2]] = Subele[,3]
                indx              = lower.tri(MAT,diag=DIAG)
                ele               = cbind(which(indx,arr.ind=T),MAT[indx])
                Nele              = dim(ele)[1]
                Subele            = ele[sample(Nele,floor(Nele*P),replace=F),]
                Subele[,3]        = runif(dim(Subele)[1],0,1)
                MAT[Subele[,1:2]] = Subele[,3]
                X                 = MAT
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            }
        }
        
    }

    return(X)
    
}

shuffle4 <- function(MAT,P=1){
    ## randomly select P% elements in MAT and then set values randonly between [0,1]
    
    X = NULL

    if( !is.null(MAT) ){
        NR = dim(MAT)[1]
        NC = dim(MAT)[2]

        if( P == 0 ){ return(MAT); }
        
        if( P == 1 ){

            if( NR == 1 ){
                X        = matrix(runif(NC,0,1),nrow=NR,ncol=NC)
                names(X) = colnames(MAT)
            } else {
                  X           = matrix(runif(NR*NC,0,1),nrow=NR,ncol=NC)
                  colnames(X) = colnames(MAT)
                  rownames(X) = rownames(MAT)  
            }

        } else {

            if( NR == 1 ){
                X         = matrix(NA,nrow=NR,ncol=NC)
                indx      = seq(1,NC,1)
                indx      = sample(indx,floor(P*NC),replace=F)
                X[1,indx] = runif(floor(P*NC),0,1)                              
                X[1,]     = ifelse(is.na(X[1,]),MAT[1,],X[1,])
                #names(X)  = colnames(MAT)
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            } else {
                X            = matrix(NA,nrow=NR,ncol=NC)
                indx         = seq(1,(NR*NC),1)
                indx         = sample(indx,floor(P*NR*NC),replace=F)
                X[indx]      = runif(floor(P*NR*NC),0,1)                 
                indx         = which(is.na(X),arr.ind=T)
                X[indx]      = MAT[indx]
                colnames(X) = colnames(MAT)
                rownames(X) = rownames(MAT)
            }
            
        }#ifelse            
    }#null        
      
    return(X)

}

zScore <- function(x, useLOG=FALSE){

    x = as.numeric(x)
    
    if( useLOG ){
        x[x<=0]=NA
        x = log(x)
    }
  
    x = (x-mean(x,na.rm=T))/sd(x,na.rm=T)

    return(x)
}

## KL divergence KL(P || Q)
kl_divergence <- function(p,q){
    p = as.numeric(p)
    q = as.numeric(q)
    if( is.na(p) || is.na(q) ){return(NA)}
    else if(p==0 || q==0){ return(NA) }
    else { return(p*log(p/q)) }
}

## Jensen-Shannon divergence JS(P || Q)
js_divergence <- function(p,q){
    p = as.numeric(p)
    q = as.numeric(q)
    m = 0.5*(p+q)
    return( (0.5*(kl_divergence(p,m)) + 0.5*(kl_divergence(q,m))) )    
}


rankDistance <- function(X,Y){

    X = as.numeric(X)
    Y = as.numeric(Y)

    if( is.na(X) && is.na(Y) ){ return(NA); }

    if( is.na(X) ){ X=1e100; }
    if( is.na(Y) ){ Y=1e100; }
    
    return( sqrt( (X-Y)^2 ) );
    
}


Stabplot <- function(DF,YLAB="",TIT="",DIR="",subDIR="",shortTIT="",groups=FALSE){

    DF = as.data.frame(DF)
    Xmin = min(as.numeric(as.vector(DF$errRate)))
    Xmax = max(as.numeric(as.vector(DF$errRate)))

    if( groups ){

        colours <- c('firebrick2','royalblue2')
        shapes  <- c(16,17)
        
        gplot = ggplot(DF,aes(x=as.numeric(as.vector(errRate)),y=as.numeric(as.vector(stab)),colour=as.vector(group)))+           
             geom_point(aes(colour=as.vector(group), shape=as.vector(group)))+
            geom_errorbar(aes(ymin=as.numeric(as.vector(CI_lower)), ymax=as.numeric(as.vector(CI_upper))),show.legend=FALSE)+
            theme(legend.key=element_blank())+
        labs(x="error Rate",y=YLAB,title=TIT)+
        scale_x_continuous(breaks = seq(Xmin, Xmax, by = 0.1))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.0)),
            axis.title.y=element_text(face="bold",size=rel(1.0)),
            axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
            legend.title=element_text(face="bold",size=rel(1.0)),
            legend.text=element_text(face="bold",size=rel(1.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
            guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(7))),
                   size   = FALSE,
                   shape  = FALSE)+
            scale_color_manual("group",breaks=levels(factor(DF$group)),values=c(colours))
        
    } else {
    
    gplot = ggplot(DF,aes(x=as.numeric(as.vector(errRate)),y=as.numeric(as.vector(stab))))+
        geom_point()+geom_errorbar(aes(ymin=as.numeric(as.vector(CI_lower)), ymax=as.numeric(as.vector(CI_upper))))+
        labs(x="error Rate",y=YLAB,title=TIT)+
        scale_x_continuous(breaks = seq(Xmin, Xmax, by = 0.1))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.0)),
            axis.title.y=element_text(face="bold",size=rel(1.0)),
            axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
            legend.title=element_text(face="bold",size=rel(1.0)),
            legend.text=element_text(face="bold",size=rel(1.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))
    }
        
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot)
    dev.off()
    
    return(gplot)
    
}


Stabdf <- function( errRate, DF, group="" ){

    N   = length(errRate)
    OUT = matrix(NA,ncol=5,nrow=N)
    colnames(OUT) = c("errRate","stab","CI_lower","CI_upper","group")
    
    k=1
    for( e in 1:N ){
        if( errRate[e] == 0 ){
            OUT[k,1] = errRate[e]
            OUT[k,2] = 1
            OUT[k,3] = 1
            OUT[k,4] = 1
            OUT[k,5] = group
        } else {
            OUT[k,1] = errRate[e]
            OUT[k,2] = as.numeric(DF[[e]]$stab)
            OUT[k,3] = as.numeric(DF[[e]]$CIlower)
            OUT[k,4] = as.numeric(DF[[e]]$CIupper)
            OUT[k,5] = group
        }
        k=k+1
    }

    return(OUT)
}


avgFeatureProb <- function(DF){

    R  = dim(DF)[1];
    N  = dim(DF)[2];

    Pc = 0;

    tmp = rep(0,R);
    for( i in 1:R ){ tmp[i] = sum(DF[i,2:N]==1)/(N-1); }

    Pc = sum(tmp)/R;

    return(Pc);
}

featureFreq <- function(DF){

    R  = dim(DF)[1];
    N  = dim(DF)[2];

    freq = rep(0,R);
    for( i in 1:R ){ freq[i] = sum(DF[i,2:N]==1); }

    return(freq);
   
}

stabilityCI <- function(DF,Pi,Kbar,stab,alpha=0.05){

    ## Ref: https://www.jmlr.org/papers/volume18/17-514/17-514.pdf
    
    F   = dim(DF)[1];     #number of features
    N   = dim(DF)[2];
    M   = (N-1);          #number of random samples
    Kf  = Kbar/F;
    
    Norm = 1/( Kf*(1-Kf) );
    
    Thetai = rep(NA,M);

    for( i in 1:(M-1) ){
        ii    = (i+1);
        ki    = sum(DF[,ii],na.rm=T);
        term1 = 0;
        for( f in 1:F ){ term1 = term1 + DF[f,ii] * Pi[f]; }
        term1 = (1/F) * term1;
        term2 = (ki*Kf)/F;
        term3 = (stab/2) * ( (2*Kf*ki)/F - ki/F - Kf + 1 );
        Thetai[ii] = Norm * ( term1 - term2 + term3 );
    }

    ## the average value of Thetai
    Thetai_avg = 1/M * sum(Thetai,na.rm=T);

    ## estimate of the variance of stab
    Thetai_var = 4/M^2 * sum( (Thetai - Thetai_avg)^2, na.rm=T);

    ## the inverse cumulative of a standard normal distribution at 1−alpha/2
    ## the inverse cumulative of a standande normal at x is calculated using qnorm(x)
    Zstar    = qnorm( (1-alpha/2) );

    upper    = stab + Zstar * sqrt( Thetai_var );
    lower    = stab - Zstar * sqrt( Thetai_var );

    return(list(var=Thetai_var,upper=upper,lower=lower))
    
}

featureStability <- function(DF, alpha=0.05){
    ## Ref: https://www.jmlr.org/papers/volume18/17-514/17-514.pdf

    F   = dim(DF)[1];     #number of features
    N   = dim(DF)[2];
    M   = (N-1);          #number of random samples
    
    
    Fn  = as.character(DF[,1]);    #feature names
    Fi  = featureFreq(DF);         #frequency that feature i is selected
    Pi  = Fi/M;                    #probability that feature i is selected
    Kbar= sum(Fi)/M;               #average number of features selected
    Kf  = Kbar/F;                  #probability of a feature begin selected
    Efi = Kf * M;                  #expected number times a feature is selected
    Sf  = (M/(M-1)) * Pi * (1-Pi); #variance of feature i being selected

    ##negative log Binomial probability for selecting each feature 
    lBin= rep(NA,F);               
    lGad= rep(NA,F);
    for( i in 1:F ){
        if( Pi[i] > 0 ){
            lBin[i] = -lbinomial(M,Fi[i],Pi[i]);
            lGad[i] = lgrad(M,Fi[i],Pi[i]);
        }
    }
    
    ## sample feature selection stability measure
    stab= ( (1/F) * sum( Sf, na.rm=T ) ); 
    stab= stab / ( Kf*(1-Kf) );
    stab= 1 - stab;

    ## confidence interval on sample feature statbility
    ci = stabilityCI(DF=DF,Pi=Pi,Kbar=Kbar,stab=stab,alpha=alpha);
    
    ## expect population feature selection stability
    Pbar  = sum(Pi)/F;
    Estab = (1/F) * sum( Pi * (1-Pi) );
    Estab = Estab / ( Pbar*(1-Pbar) );
    Estab = 1 - Estab;    
    
    ##feature robustness measure
    Sigma  = sqrt( Kf*(1-Kf)*M );
    Zscore = rep(NA,F);
    if( Sigma != 0 ){ Zscore = (Fi - Efi)/Sigma; }
    
    CN  = c("feature","Fi","Pi","Sf","nlBin","lgrad","robst")
    OUT = matrix(NA,ncol=length(CN),nrow=F);
    colnames(OUT) = CN
    OUT     = as.data.frame(OUT)
    OUT[,1] = as.character(Fn)
    OUT[,2] = as.numeric(Fi)
    OUT[,3] = as.numeric(Pi)
    OUT[,4] = as.numeric(Sf)
    OUT[,5] = as.numeric(lBin)
    OUT[,6] = as.numeric(lGad)
    OUT[,7] = as.numeric(Zscore)

    return(list(stab=stab,Estab=Estab, Vstab=ci$var, CIupper=ci$upper, CIlower=ci$lower,
                alpha=alpha, OUT=OUT));
    
}

#fullSetRobustness <- function(){}

featureRobustness <- function(DF){
# Ref: http://memetic-computing.org/publication/journal/MBGA.pdf
    
    F   = dim(DF)[1];     #number of features
    N   = dim(DF)[2];
    M   = (N-1);          #number of random samples
    
    Fi  = as.character(DF[,1])          #feature names
    Si  = featureFreq(DF) #frequency that feature i is selected
    Sc  = sum( (Si/M) );  #average number of features selected
    
    Psi = Sc/F;           #probability of a feature begin selected

    Esi = Psi * M;        #expected number of times a feature is selected

    Sigma = sqrt( Psi*(1-Psi)*M ); #the standard deviation

    CN = c("feature","robustness")
    Zscore = matrix(NA,ncol=length(CN),nrow=F);
    colnames(Zscore) = CN
    Zscore[,1]       = as.character(Fi)

    if( Sigma != 0 ){ Zscore[,2] = (Si - Esi) / Sigma; }

    return(Zscore);
    
}


linearRankSel <- function(X,decreasing=FALSE, negative=FALSE){
# Ref: https://stackoverflow.com/questions/20290831/how-to-perform-rank-based-selection-in-a-genetic-algorithm

    F = dim(X)[1]

    if( !negative ){
        Fsub = X[as.numeric(X[,2]) > 0,]
    } else {
        Fsub = X
    }
        
    indx = order(as.numeric(Fsub[,2]),decreasing=!decreasing)

    Pndx = order(as.numeric(Fsub[,2]),decreasing=decreasing)
    
    R    = dim(Fsub)[1]

    Rank = seq(1,R,1)
    
    SumR = (R+1)*(R/2)

    Pf   = Rank/SumR

    CN   = c(colnames(X),"rank","prob")
    out  = matrix(NA,ncol=length(CN),nrow=F)
    colnames(out) = CN
    out[,1]       = as.character(X[,1])
    out[,2]       = as.character(X[,2])
    
    out[match(Fsub[indx,1],out[,1]),3] = Rank
    out[match(Fsub[Pndx,1],out[,1]),4] = Pf

    return(out)
    
}

buildPhenotypeProbMatrix <- function(VAR, ii, CT){

    VALS = names(table(VAR))
    xx   = CT[match(VAR,VALS),ii]
    return(xx/sum(xx,na.rm=T))
    
    }

buildCommunitySizeTable <- function(VATTS, ID="SPARKID", ALG=NULL ){

    RES = NULL
        
    if( !is.null(ALG) ){

        IDindx = which(colnames(VATTS)==ID)
        Ccindx = which(colnames(VATTS)==ALG)

        if( length(IDindx) !=0  && length(Ccindx) != 0 ){
        
            Cmax   = max(VATTS[,Ccindx])
            C      = seq_along(1:Cmax)

            RES    = matrix(0,Cmax,2)
            RES[,1]= C
            for( c in 1:Cmax ){
                RES[c,2] = sum(VATTS[,Ccindx]==c)
            }
            
        }
    }

    return(RES)
    
}

buildPhenotypeMatrix <- function( VATTS, FIN, ID="SPARKID" ){

    RES = NULL
    FN  = length(names(FIN))
    
    if( !is.null(ID) ){

        IDindx = which(colnames(VATTS)==ID)

        if( length(IDindx) !=0  ){           

            IDS <- as.vector(VATTS[,IDindx])
            RN  <- length(IDS)
            CN  <- c()

            for( f in 1:FN ){
                N   = length(colnames(FIN[[f]]))
                ii  = seq(3,N,1)
                CN <- c(CN,sprintf("%s_%s",names(FIN)[f],colnames(FIN[[f]])[ii]))
            }

            RES <- matrix(NA,RN,length(CN))
            colnames(RES) <- CN

            for( f in 1:FN ){

                N   = length(colnames(FIN[[f]]))
                ii  = seq(3,N,1)
                CN  = sprintf("%s_%s",names(FIN)[f],colnames(FIN[[f]])[ii])
                id  = as.character(unlist(FIN[[f]][,.SD,.SDcol=1]))

                Rindx = match(IDS,id)

                for( i in 1:length(ii) ){
                    val   = as.vector(unlist(FIN[[f]][Rindx,.SD,.SDcol=ii[i]]))                    
                    Cindx = match(CN[i],colnames(RES))
                    RES[,Cindx] = val
                }
                
            }

        }
    }

    return(RES)

}
            
buildPhenotypeCorrelationMatrix <- function( MAT ){

    COR  = NULL
    CORz = NULL

    if( !is.null(MAT) ){
    
    CN   = colnames(MAT)
     N   = length(CN)
    COR  = matrix(NA,N,N)
    CORz = matrix(NA,N,N)

    colnames(COR)  = CN
    rownames(COR)  = CN
    colnames(CORz) = CN
    rownames(CORz) = CN
    
            for( i in 1:N ){
                for( j in 1:N ){

                    if( i <= j ){
                        
                        x    = MAT[,i]
                        y    = MAT[,j]

                        indx = !is.na(x) & !is.na(y)

                        val  = cor(x[indx],y[indx])
                        COR[i,j] = val
                        COR[j,i] = val

                        xz   = x[indx]
                        yz   = y[indx]

                        xz   = (xz - mean(xz))/sd(xz)
                        yz   = (yz - mean(yz))/sd(yz)
                        
                        val  = cor(xz,yz)
                        CORz[i,j] = val
                        CORz[j,i] = val
                        
                    }                    
                }
            }
            
        }
    
    
     return(list(COR=COR,CORz=CORz))
    
}

buildPhenotypeNMImatrix <- function( MAT ){

    NMI  = NULL

    if( !is.null(MAT) ){
    
    CN  = colnames(MAT)
     N  = length(CN)
    NMI  = matrix(NA,N,N)


    colnames(NMI)  = CN
    rownames(NMI)  = CN
    
            for( i in 1:N ){
                for( j in 1:N ){

                    if( i <= j ){                        
                   
                        res  = MutualInformation(MAT[,i],MAT[,j])
                        NMI[i,j] = res$NMI
                        NMI[j,i] = res$NMI
                        
                    }                    
                }
            }
       
    }   
    
    return(list(NMI=NMI))
    
}

extractSubMatrix <- function( MAT, TARGETS ){

    RES <- NULL
    
    if( !is.null(MAT) && !is.null(TARGETS) ){
        indx = match(TARGETS, colnames(MAT))
        if( length(indx) != 0 ){
            RES = MAT[indx,indx]
            #colnames(RES) = colnames(MAT)
            #rownames(RES) = colnames(MAT)
        }
    }

    return(RES)
    
}

saveNMImatrix <- function( MAT, DIR, subDIR, TITLE="" ){

    cat("Save ", TITLE,"...")
    
    saveRDS(MAT, sprintf("%s/%s/%s.RDS",DIR,subDIR,TITLE),compress=TRUE)

    cat("done.\n")
}


saveMImatrices <- function( MAT, GROUPS, DIR, subDIR, TITLE="" ){

    cat("building Conditional MI ", TITLE,"...")
    
    ## build I(X,Y|Z) matrix for Z = GROUPS
    res = buildPhenotypeIImatrix(MAT,GROUPS)

    cat("done! ")
    
    saveRDS(res, sprintf("%s/%s/%s.RDS",DIR,subDIR,TITLE),compress=TRUE)

    cat("Saved ", TITLE, ".\n")
}

CSFplot <- function(DF,YLAB="",TIT="",DIR="",subDIR="",shortTIT=""){

    DF = as.data.frame(DF)
    Xmin = min(as.numeric(as.vector(DF$errRate)))
    Xmax = max(as.numeric(as.vector(DF$errRate)))
    gplot = ggplot(DF,aes(x=as.numeric(as.vector(errRate)),y=as.numeric(as.vector(median))))+
        geom_point()+geom_errorbar(aes(ymin=as.numeric(as.vector(min)), ymax=as.numeric(as.vector(max))))+
        labs(x="error Rate",y=YLAB,title=TIT)+
        scale_x_continuous(breaks = seq(Xmin, Xmax, by = 0.1))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.0)),
            axis.title.y=element_text(face="bold",size=rel(1.0)),
            axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
            legend.title=element_text(face="bold",size=rel(1.0)),
            legend.text=element_text(face="bold",size=rel(1.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))

    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot)
    dev.off()
    
    return(gplot)
    
}

CFSdf <- function( errRate, DF, INDX=1 ){

    N   = length(errRate)
    OUT = matrix(NA,ncol=5,nrow=N)
    colnames(OUT) = c("errRate","median","min","max","max.ind")
    
    k=1
    for( e in 1:N ){
        if( errRate[e] == 0 ){
            OUT[k,1] = errRate[e]
            OUT[k,2] = DF[[e]][INDX]
            OUT[k,3] = DF[[e]][INDX]
            OUT[k,4] = DF[[e]][INDX]
            OUT[k,5] = 2
        } else {
            OUT[k,1] = errRate[e]
            OUT[k,2] = median(DF[[e]][,INDX])
            OUT[k,3] = min(DF[[e]][,INDX])
            OUT[k,4] = max(DF[[e]][,INDX])
            OUT[k,5] = (which.max(DF[[e]][,INDX])+1)
        }
        k=k+1
    }

    return(OUT)
}

CFS <- function(X=NULL, fNMI=NULL, cNMI=NULL ){
    ## Correlation feature selection, using symmetric uncertainty (i.e. NMI). 
    ## Ref: https://en.wikipedia.org/wiki/Feature_selection
    ## X    ==> feature selection vector,
    ## X[i] == 1 if feature selected and X[i] == 0 if feature not selected
    ## fNMI ==> feature-feature correlation
    ## cNMI ==> feature-class   correlation

    CFS = NA
    
    if( !is.null(X) && !is.null(fNMI) && !is.null(cNMI) ){

        X[is.na(X)]=0
        cNMI[1,is.na(cNMI[1,])]=0
        fNMI[is.na(fNMI)]=0
        diag(fNMI)=0
        
        num   = cNMI[1,] %*% X;

        term1 = sum( X );
        term2 = 2 * X %*% fNMI %*% X; 

        CFS = num/(term1 + term2);
        
    }


    return(CFS)
}

CCFS <- function(X=NULL, fNMI=NULL, cNMI=NULL, CNMI=NULL ){
    ## Conditional Correlation feature selection, using symmetric uncertainty (i.e. NMI). 
    ## Ref: https://en.wikipedia.org/wiki/Feature_selection
    ## X    ==> feature selection vector,
    ## X[i] == 1 if feature selected and X[i] == 0 if feature not selected
    ## fNMI ==> feature-feature correlation
    ## cNMI ==> feature-class   correlation
    ## CNMI ==> feature-feature given class correlation
    
    CCFS = NA
    
    if( !is.null(X) && !is.null(fNMI) && !is.null(cNMI) && !is.null(CNMI) ){

        X[is.na(X)]=0
        cNMI[1,is.na(cNMI[1,])]=0
        fNMI[is.na(fNMI)]=0
        diag(fNMI)=0

        CNMI[is.na(CNMI)]=0
        diag(CNMI)=0
        
        num   = cNMI[1,] %*% X;

        term1 = sum( X );
        term2 = 2 * X %*% (fNMI-CNMI) %*% X; 

        CCFS = num/(term1 + term2);
        
    }


    return(CCFS)
}


getNextFi <- function( Fsel, Fi ){

    indx = which(Fsel==Fi)

    if( length(indx) == 0 ){ return(NULL); }
    else { return(Fsel[indx+1]); }
    
}
    
addFi <- function( Fset, Fi ){

    if( is.null(Fset) ){ return(Fi); }
    
    indx = which(Fset==Fi)

    if( length(indx) == 0 ){ return(c(Fi,Fset)); }
    else                   { return(Fset);       }
}

removeFi <- function( Fsel, Fi ){

    indx = which(Fsel==Fi)
    if( length(indx) == 0 ){ return(Fsel) }
    else { return(Fsel[-indx]); }
}
    
FCCFrule <- function( Fi, Fj, fSU, cSU, pSU, CSU, NCLASSES ){
    ## Fast Class Correlation Filter rule (FCCF)       

    ## Does the feature Fj form a Markov Blacket around feature Fi (in other words,
    ## is feature Fi subsume/absorb by feature Fj), if so
    ## feature Fi is redundat and we can remove it.

    ## return 1 => remove feature
    ## return 0 => don't remove feature
    
    ii = which(colnames(pSU)==Fi)
    jj = which(colnames(pSU)==Fj)

    if( length(ii) == 0 || length(jj) == 0 ){ return(0); }   

    if( is.na(fSU[ii,jj]) || is.na(cSU[ii]) ){ return(0); }
    
    perClass = pSU[,jj] >= pSU[,ii] 
    
    if( (sum(perClass,na.rm=T) == NCLASSES) && (fSU[ii,jj] >= cSU[ii]) ){
        ## feature Fi redundant and can be removed.
        return(1);
    }
    
    ## feature Fi not redundant so keep.
    return(0);
    
}

FCBFrule <- function( Fi, Fj, fSU, cSU, pSU, CSU, NCLASSES ){
    ## Fast Correlation Based Filter rule (FCCF)       

    ## Does the feature Fj form a Markov Blacket around feature Fi (in other words,
    ## is feature Fi subsume/absorb by feature Fj), if so
    ## feature Fi is redundat and we can remove it.

    ## return 1 => remove feature
    ## return 0 => don't remove feature

    #cat("Fi: ", Fi, " Fj: ", Fj, "\n")
    
    ii = which(colnames(pSU)==Fi)
    jj = which(colnames(pSU)==Fj)

    #cat("ii: ", ii, " jj: ", jj, "\n")
    #cat("cSU[ii]: ", cSU[ii], " fSU: ", fSU[ii,jj], " cSU[jj]: ", cSU[jj], "\n")
    
    if( length(ii) == 0 || length(jj) == 0 ){ return(0); }   

    if( is.na(cSU[ii]) || is.na(cSU[jj]) || is.na(fSU[ii,jj]) ){ return(0); }
    
    if( (cSU[jj] >= cSU[ii]) && (fSU[ii,jj] >= cSU[ii]) ){
        ## feature Fi redundant and can be removed.
        return(1);
    }

    ## feature Fi not redundant so keep.
    return(0);
    
}
 

Jcmirule <- function( Fi, Fj, fSU, cSU, pSU, CSU, NCLASSES ){
    ## Jcmi rule

    ## Does the feature Fj form a Markov Blacket around feature Fi (in other words,
    ## is feature Fi subsume/absorb by feature Fj), if so
    ## feature Fi is redundat and we can remove it.

    ## Jcmi(Fi) = NMI(Xi;Y) - NMI(Fi;Fj) + NMI(Fi;Fj|Y)
    ## NMI(Fi;Fj) - NMI(Fi;Fj|Y) >= NMI(Xi;Y)
    
    ## return 1 => remove feature
    ## return 0 => don't remove feature
    
    ii = which(colnames(pSU)==Fi)
    jj = which(colnames(pSU)==Fj)

    if( length(ii) == 0 || length(jj) == 0 ){ return(0); }
    
    ## expected feature correlation condition on class
    #condSU = 0
    #for( i in 1:NCLASSES ){ condSU = condSU + CSU[[i]][ii,jj]; }
    #condSU = condSU/NCLASSES;
    condSU = CSU[ii,jj]
    ## -----

    if( is.na(cSU[ii]) || is.na(condSU) || is.na(fSU[ii,jj]) ){ return(0); }

    ## let's keep this... matches condition in original FCBF rule
    if( (cSU[jj] >= cSU[ii]) && (fSU[ii,jj] - condSU) >= cSU[ii] ){

    ## old condition, which is also fine.
    ##if( (fSU[ii,jj] - condSU) >= cSU[ii] ){
    
        ## feature Fi redundant and can be removed.
        return(1);
    }

    ## feature Fi not redundant so keep.
    return(0);
    
}

CFCCFrule <- function( Fi, Fj, fSU, cSU, pSU, CSU, NCLASSES ){
    ## Conditional Fast Class Correlation Filter rule (CFCCF)       

    ## Does the feature Fj form a Markov Blacket around feature Fi (in other words,
    ## is feature Fi subsume/absorb by feature Fj), if so
    ## feature Fi is redundat and we can remove it.

    ## return 1 => remove feature
    ## return 0 => don't remove feature
    
    ii = which(colnames(pSU)==Fi)
    jj = which(colnames(pSU)==Fj)

    if( length(ii) == 0 || length(jj) == 0 ){ return(0); }

    
    perClass = pSU[,jj] >= pSU[,ii] 

    #condSU = 0
    #for( i in 1:NCLASSES ){ condSU = condSU + CSU[[i]][ii,jj]; }
    #condSU = condSU/NCLASSES;    
    condSU = CSU[ii,jj]
    
    if( (sum(perClass) == NCLASSES) && (fSU[ii,jj] >= (cSU[ii] + condSU)) ){
        ## feature Fi redundant and can be removed.
        return(1);
    }

    ## feature Fi not redundant so keep.
    return(0);
    
}


##REFS:
## [1]: FEATURE SELECTION VIA JOINT LIKELIHOOD, 2012 (http://www.cs.man.ac.uk/~gbrown/publications/pocockPhDthesis.pdf).
## [2]: Efficient Feature Selection via Analysis ofRelevance and Redundancy. Lei Yu & Huan Liu, Journal of Machine Learning Research 5 (2004) 1205–1224 (https://jmlr.org/papers/volume5/yu04a/yu04a.pdf)
## [3]: Scalable Feature Selection for Multi-class Problems, Boris Chidlovskii and Loıc Lecerf, 2008 (https://link.springer.com/content/pdf/10.1007%2F978-3-540-87479-9_33.pdf)

FCBF <- function( fNMI, cNMI, pNMI, CNMI, delta=NULL, excSET=NULL, print=FALSE, rule="FCCF", runSel=FALSE ){
## Fast Correlation-Based Filter Algorithm (FCBF) 
    
    Fout <- NULL
        
 ## Fast Class Correlation Filter
 ## fNMI  => feature-by-feature NMI
 ## cNMI  => class-by-feature NMI
 ## pNMI  => per class-by-feature NMI
 ## delta => irrelevance threshold 
    if( !is.null(fNMI) && !is.null(cNMI) && !is.null(pNMI) ){

        #if( (rule == 3 || rule == 4) && is.null(CNMI) ){ break; }

        if( print ){ cat("using rule", rule, "\n"); } 

        if( runSel == TRUE ){
            if( (rule == "FCBF") || (rule == "FCCF") ){
                score = variableSelectionFCBF( cNMI=cNMI, fNMI=fNMI );
            }
        
            if( (rule == "Jcmi") || (rule == "CFCCF") ){
                score = variableSelectionJcmi( cNMI=cNMI, fNMI=fNMI, CNMI=CNMI );
            }
        }
        
        CLASSES = dim(pNMI)[1] ## number of classes
        N       = dim(pNMI)[2] ## number of features
        Irr     = c()
        Rel     = seq_along(1:N)
        
        sig     = data.frame(a=as.character(),b=as.numeric());
        
        ## step 1, select subset of relevant features
        if( is.null(excSET) ){

            
        ## select subset of relevant features using cut-off on SU(Y,Fi) values.
        if( is.null(delta) ){ delta = quantile(pNMI)[2]; }
        
        RELindx = cNMI > delta

        FpNMI = pNMI[,RELindx]
        if( is.vector(cNMI) ){ FcNMI = cNMI[RELindx]; }
        if( is.matrix(cNMI) ){ FcNMI = cNMI[,RELindx]; }
        FfNMI = fNMI[RELindx,RELindx]

        if( !is.null(CNMI) ){ CNMI = CNMI[RELindx,RELindx]; }            
        #if( !is.null(CNMI) ){
        #    for( i in 1:length(names(CNMI)) ){ CNMI[[i]] = CNMI[[i]][RELindx,RELindx]; }
        #}
            
        tmp   = data.frame( a=as.character(colnames(pNMI)[!RELindx]),
                            b=as.numeric(rep(-1,sum(!RELindx))) );
        sig   = rbind(sig,tmp);

        Irr   = !RELindx
        Rel   = RELindx    
            
        } else {

            ## user defined sub set of features to exclude
            EXCindx = match(excSET,colnames(pNMI))
            EXCindx = EXCindx[!is.na(EXCindx)]

            if( length(EXCindx) != 0 ){
            
                FpNMI = pNMI[,-EXCindx]
                if( is.vector(cNMI) ){ FcNMI = cNMI[-EXCindx]; }
                if( is.matrix(cNMI) ){ FcNMI = cNMI[,-EXCindx]; }
                FfNMI = fNMI[-EXCindx,-EXCindx]

                if( !is.null(CNMI) ){ CNMI = CNMI[-EXCindx,-EXCindx]; }
                #if( !is.null(CNMI) ){
                #    for( i in 1:length(names(CNMI)) ){ CNMI[[i]] = CNMI[[i]][-EXCindx,-EXCindx]; }
                #}              
                
                Irr = EXCindx
                Rel = Rel[-EXCindx]
                
            } else {
                FpNMI = pNMI
                FcNMI = cNMI
                FfNMI = fNMI
            }
            
            tmp   = data.frame( a=as.character(colnames(pNMI)[EXCindx]),
                                b=as.numeric(rep(-1,length(EXCindx))) );
            sig   = rbind(sig,tmp);
            
        }

    
        ## step 2, select predominat features from relevant ones,
        ## i.e. remove redundant features which have a Markov Blanket.

        if( is.vector(FcNMI)){ Fsel = names(FcNMI);    }
        if( is.matrix(FcNMI)){ Fsel = colnames(FcNMI); }
        
        ## order features
        Frank = order(FcNMI,decreasing=T)
        Fsel  = Fsel[Frank]
        Frank = cbind(Fsel,seq(1,length(Fsel),1))       
        
        ## get first feature
        Fj = Fsel[1]

        while( !is.null(Fj) ){

            ##Fopt  = addFi(Fpivot, Fopt);
            Fi = getNextFi(Fsel=Fsel, Fi=Fj);

            while( !is.null(Fi) ){

                switch(rule,
                       "FCBF"={
                  signal = FCBFrule(Fi=Fi, Fj=Fj, fSU=FfNMI, cSU=FcNMI, pSU=FpNMI, CSU=CNMI, NCLASSES=CLASSES);
                       },
                       "FCCF"={
                  signal = FCCFrule(Fi=Fi, Fj=Fj, fSU=FfNMI, cSU=FcNMI, pSU=FpNMI, CSU=CNMI, NCLASSES=CLASSES);
                       },
                       "Jcmi"={
                  signal = Jcmirule(Fi=Fi, Fj=Fj, fSU=FfNMI, cSU=FcNMI, pSU=FpNMI, CSU=CNMI, NCLASSES=CLASSES);
                       },
                       "CFCCF"={
                  signal = CFCCFrule(Fi=Fi, Fj=Fj, fSU=FfNMI, cSU=FcNMI, pSU=FpNMI, CSU=CNMI, NCLASSES=CLASSES);          
                       })
                
                if( signal == 1 ){                    
                    Fsel = removeFi(Fsel=Fsel, Fi=Fi);                
                    if( print ){ cat("remove: ", Fi, "\n"); }
                    tmp = data.frame(a=as.character(Fi),b=as.numeric(signal) );
                    sig = rbind(sig,tmp);
                }                
                Fi = getNextFi(Fsel=Fsel, Fi=Fi);
            }
            Fj = getNextFi(Fsel=Fsel, Fi=Fj);
        }

    tmp = data.frame(a=as.character(Fsel),b=as.numeric(rep(2,length(Fsel))));
    sig = rbind(sig,tmp);
        

    ##output 
    CN = c("feature","irrelevant","relevant","predominat","rank","score","signal")
    Fout = matrix(0,ncol=length(CN),nrow=N)
    colnames(Fout)                  = CN
    Fout[,1]                        = colnames(pNMI)    
    if( length(Irr) != 0 ){ Fout[Irr,2] = 1; }
    if( length(Rel) != 0 ){ Fout[Rel,3] = 1; }
    #Fout[!RELindx,2]                = 1
    #Fout[RELindx,3]                 = 1
    Fout[match(Fsel,Fout[,1]),4]        = 1;
    if( runSel ){     
        Fout[match(score[,1],Fout[,1]),5]   = score[,4];
        Fout[match(score[,1],Fout[,1]),6]   = score[,5];
    } else {
        Fout[,5]                            = rep(NA,length(Fout[,5]));
        Fout[,6]                            = rep(NA,length(Fout[,6]));
    }
    #Fout[match(Frank[,1],Fout[,1]),5]   = Frank[,2];
    #Fout[,5] = ifelse( Fout[,6] == 0, NA, Fout[,6]);
    Fout[match(sig[,1],Fout[,1]),7]     = sig[,2];
    }
        
    return(Fout);
    
}
    
symmetricUncertainty <- function( MAT, GROUPS ){

    NMI  = NULL
    cNMI = NULL
    
    if( !is.null(MAT) ){
    
    CN   = colnames(MAT)
    RN   = colnames(GROUPS)    
    Ncol = length(CN)
    Nrow = length(RN)
    NMI  = matrix(NA,Nrow,Ncol)

    colnames(NMI)  = CN
    rownames(NMI)  = RN
    
            for( i in 1:Nrow ){
                for( j in 1:Ncol ){

                        res  = MutualInformation(GROUPS[,i],MAT[,j])
                        NMI[i,j] = res$NMI

                
            }
            
        }

        cNMI = matrix(NA,1,Ncol)
        colnames(cNMI) = CN
        rownames(cNMI) = "ALL"

        if( Nrow == 1 ){
            cNMI           = NMI
            rownames(cNMI) = "ALL"
        } else {
        
            for( i in 1:(Nrow-1) ){
                GROUPS[,(i+1)] = GROUPS[,(i+1)] + max(GROUPS[,i])
            }

            a = GROUPS[,1]
            for( i in 2:Nrow ){
                a = c(a,GROUPS[,i])
            }

            for( j in 1:Ncol ){        
                b=rep(MAT[,j],Nrow) 
                res=MutualInformation(a,b)
                cNMI[1,j] = res$NMI
            }

        }#ifelse
    }
    
     return(list(NMI=NMI, cNMI=cNMI))
    
    
}

variableSelection <- function( VARS, Y, RES ){

    S  = NULL
    Sn = NULL
    
    ## step 1: calculate I(X_k;Y) for each variable

    temp  = list()
    Nf    = length(colnames(VARS))
    F     = seq_along(1:Nf)

    S     = matrix(NA,ncol=5,nrow=Nf)
    S[,1] = colnames(VARS)
    S[,2] = rep(FALSE,Nf)
    S[,3] = rep(TRUE,Nf)
    colnames(S) <- c("variable name","S","!S","rank","Jcmi")

    Sn     = matrix(NA,ncol=5,nrow=Nf)
    Sn[,1] = colnames(VARS)
    Sn[,2] = rep(FALSE,Nf)
    Sn[,3] = rep(TRUE,Nf)
    colnames(Sn) <- c("variable name","S","!S","rank","Jcmi")
    
    for( k in 1:Nf ){
        temp[[k]]      <- MutualInformation(VARS[,k],Y)
        names(temp)[k] <- colnames(VARS)[k]
    }

    MIXkY            <- cbind( sapply(temp,Extract,ii=4) ) 
    NMIXkY           <- cbind( sapply(temp,Extract,ii=5) ) 
    rownames(MIXkY)  <- NULL
    rownames(NMIXkY) <- NULL

    ## step 2: select variable with highest I(X_k;Y) and put into set S
    Fadd1      = which.max(as.numeric(MIXkY))
    S[Fadd1,2] = TRUE
    S[Fadd1,3] = FALSE
    S[Fadd1,4] = 1
    S[Fadd1,5] = round(max(as.numeric(MIXkY),na.rm=T),3) #round(res$MI[Fadd1,Fadd1],3)

    Fnadd1       = which.max(as.numeric(NMIXkY))
    Sn[Fnadd1,2] = TRUE
    Sn[Fnadd1,3] = FALSE
    Sn[Fnadd1,4] = 1
    Sn[Fnadd1,5] = round(max(as.numeric(NMIXkY),na.rm=T),3) #round(res$MI[Fadd1,Fadd1],3)
    
    
    ## step 3: keep selecting the variable with the next highest Jcmi score, and add to set S 
    for( i in 1:(Nf-1) ){

        jcmi = Jcmi(S=S, IXkY=MIXkY, IXjXk=RES$IXY, IXjXkCY=RES$MI)

        Jmax = max(jcmi,na.rm=T)    
        Fadd = which.max(as.numeric(jcmi))
        if( !is.na(jcmi[Fadd]) && length(Fadd) > 0 ){
            S[Fadd,2] = TRUE
            S[Fadd,3] = FALSE
            S[Fadd,4] = (i+1)
            S[Fadd,5] = round(Jmax,3)
        }

        jcmi = Jcmi(S=Sn, IXkY=NMIXkY, IXjXk=RES$NMIXY, IXjXkCY=RES$NMI)

        Jmax = max(jcmi,na.rm=T)    
        Fadd = which.max(as.numeric(jcmi))
        if( !is.na(jcmi[Fadd]) && length(Fadd) > 0 ){
            Sn[Fadd,2] = TRUE
            Sn[Fadd,3] = FALSE
            Sn[Fadd,4] = (i+1)
            Sn[Fadd,5] = round(Jmax,3)
        }

        
    }

    return(list(S=S,Sn=Sn))
    
}

variableSelectionJcmi <- function( cNMI, fNMI, CNMI ){

    Sn = NULL
    
    ## step 1: calculate I(X_k;Y) for each variable
    temp  = list()
    Nf    = length(colnames(cNMI))
    F     = seq_along(1:Nf)   

    Sn     = matrix(NA,ncol=5,nrow=Nf)
    Sn[,1] = colnames(cNMI)
    Sn[,2] = rep(FALSE,Nf)
    Sn[,3] = rep(TRUE,Nf)
    colnames(Sn) <- c("variable name","S","!S","rank","Jcmi")
        

    ## step 2: select variable with highest I(X_k;Y) and put into set Sn  
    Fnadd1       = which.max(as.numeric(cNMI))
    Sn[Fnadd1,2] = TRUE
    Sn[Fnadd1,3] = FALSE
    Sn[Fnadd1,4] = 1
    Sn[Fnadd1,5] = round(max(as.numeric(cNMI),na.rm=T),3) #round(res$MI[Fadd1,Fadd1],3)
    
    
    ## step 3: keep selecting the variable with the next highest Jcmi score, and add to set S 
    for( i in 1:(Nf-1) ){

        jcmi = Jcmi(S=Sn, IXkY=cNMI, IXjXk=fNMI, IXjXkCY=CNMI)       

        Jmax = max(jcmi,na.rm=T)    
        Fadd = which.max(as.numeric(jcmi))
        if( !is.na(jcmi[Fadd]) && length(Fadd) > 0 ){
            Sn[Fadd,2] = TRUE
            Sn[Fadd,3] = FALSE
            Sn[Fadd,4] = (i+1)
            Sn[Fadd,5] = round(Jmax,3)
        }

        
    }

    return(Sn=Sn)
    
}

variableSelectionFCBF <- function( cNMI, fNMI ){

    Sn = NULL
    
    ## step 1: calculate I(X_k;Y) for each variable
    temp  = list()
    Nf    = length(colnames(cNMI))
    F     = seq_along(1:Nf)   

    Sn     = matrix(NA,ncol=5,nrow=Nf)
    Sn[,1] = colnames(cNMI)
    Sn[,2] = rep(FALSE,Nf)
    Sn[,3] = rep(TRUE,Nf)
    colnames(Sn) <- c("variable name","S","!S","rank","fcbf")
        

    ## step 2: select variable with highest I(X_k;Y) and put into set Sn  
    Fnadd1       = which.max(as.numeric(cNMI))
    Sn[Fnadd1,2] = TRUE
    Sn[Fnadd1,3] = FALSE
    Sn[Fnadd1,4] = 1
    Sn[Fnadd1,5] = round(max(as.numeric(cNMI),na.rm=T),3) #round(res$MI[Fadd1,Fadd1],3)
    
    
    ## step 3: keep selecting the variable with the next highest Jcmi score, and add to set S 
    for( i in 1:(Nf-1) ){

        fcbf = Fcbf(S=Sn, IXkY=cNMI, IXjXk=fNMI)       

        Jmax = max(fcbf,na.rm=T)    
        Fadd = which.max(as.numeric(fcbf))
        if( !is.na(fcbf[Fadd]) && length(Fadd) > 0 ){
            Sn[Fadd,2] = TRUE
            Sn[Fadd,3] = FALSE
            Sn[Fadd,4] = (i+1)
            Sn[Fadd,5] = round(Jmax,3)
        }

        
    }

    return(Sn=Sn)
    
}


buildPhenotypeIImatrix <- function( MAT, GROUPS ){

    MI   =NULL
    II   =NULL
    NMI  =NULL
    IXY  =NULL
    IXZ  =NULL
    NMIXY=NULL
    NMIXZ=NULL
    
    if( !is.null(MAT) ){
    
        CN    = colnames(MAT)
        N     = length(CN)
        MI    = matrix(NA,N,N)
        II    = matrix(NA,N,N)
        NMI   = matrix(NA,N,N)
        IXY   = matrix(NA,N,N)
        IXZ   = matrix(NA,N,N)
        NMIXY = matrix(NA,N,N)
        NMIXZ = matrix(NA,N,N)
        
        colnames(MI)  = CN
        rownames(MI)  = CN

        colnames(II)  = CN
        rownames(II)  = CN

        colnames(NMI) = CN
        rownames(NMI) = CN

        colnames(IXY) = CN
        rownames(IXY) = CN

        colnames(IXZ) = CN
        rownames(IXZ) = CN

        colnames(NMIXY) = CN
        rownames(NMIXY) = CN

        colnames(NMIXZ) = CN
        rownames(NMIXZ) = CN
        
        for( i in 1:N ){
            for( j in 1:N ){

                if( i <= j ){                        
                   
                    res  = conditionalMI(MAT[,i],MAT[,j], GROUPS)

                    MI[i,j]  = res$MIabCc
                    MI[j,i]  = res$MIabCc

                    II[i,j]  = res$IIabc
                    II[j,i]  = res$IIabc

                    NMI[i,j] = res$NMI
                    NMI[j,i] = res$NMI

                    IXY[i,j] = res$MIxy$MI
                    IXY[j,i] = res$MIxy$MI

                    IXZ[i,j] = res$MIxz$MI
                    IXZ[j,i] = res$MIxz$MI

                    NMIXY[i,j] = res$MIxy$NMI
                    NMIXY[j,i] = res$MIxy$NMI

                    NMIXZ[i,j] = res$MIxz$NMI
                    NMIXZ[j,i] = res$MIxz$NMI
                    
                }                    
            }
        }
        
    }
    
    
    return(list(MI=MI, NMI=NMI, II=II, IXY=IXY, NMIXY=NMIXY, IXZ=IXZ, NMIXZ=NMIXZ))
    
}


buildCommunityPhenotypeTable <- function( VATTS, FIN, ID="SPARKID", ALG=NULL, normROW=TRUE ){

    RES = NULL
    RAW = NULL
    FN  = length(names(FIN))
    
    if( !is.null(ALG) ){

        IDindx = which(colnames(VATTS)==ID)
        Ccindx = which(colnames(VATTS)==ALG)

        if( length(IDindx) !=0  && length(Ccindx) != 0 ){
        
            Cmax   = max(VATTS[,Ccindx])
            C      = seq_along(1:Cmax)

            OO  <- list()
            RW  <- list()
            TOT <- list()

            for( f in 1:FN ){

                N   = length(colnames(FIN[[f]]))
                ii  = seq(3,N,1)

                oo           <- matrix(0,Cmax,length(ii))
                colnames(oo) <- sprintf("%s_%s",names(FIN)[f],colnames(FIN[[f]])[ii])

                tot          <- matrix(0,1,length(ii))
                colnames(tot)<- sprintf("%s_%s",names(FIN)[f],colnames(FIN[[f]])[ii])

                for( i in 1:length(ii) ){
                    val      = as.vector(unlist(FIN[[f]][,.SD,.SDcol=ii[i]])) 
                    tot[1,i] = sum(val,na.rm=T)
                }
    
                OO[[f]]   <- oo
                names(OO) <- names(FIN)[f]

                RW[[f]]   <- oo
                names(RW) <- names(FIN)[f]
                
                TOT[[f]]  <- tot
                names(TOT)<- names(FIN)[f]

            }


            for( f in 1:FN ){

                N  = length(colnames(FIN[[f]]))
                ii = seq(3,N,1)
    
                for( c in 1:Cmax ){

                    indx   = ifelse(!is.na(match(VATTS[,Ccindx],C[c])),TRUE,FALSE)
                    cids   = as.character(VATTS[indx,IDindx])

                    id    =  as.character(unlist(FIN[[f]][,.SD,.SDcol=1]))
                    Findx = ifelse(!is.na(match(id,cids)),TRUE,FALSE)

                    for( i in 1:length(ii) ){    
                        val             = as.vector(unlist(FIN[[f]][Findx,.SD,.SDcol=ii[i]]))
                        Raw             = sum(val != 0 & !is.na(val))
                        Sum             = sum(val,na.rm=T)
                        Norm            = TOT[[f]][1,i]
                        term            = 0
                        if( Norm != 0 ){ term = Sum/Norm }

                        RW[[f]][C[c],i] = RW[[f]][C[c],i] + Raw
                        OO[[f]][C[c],i] = OO[[f]][C[c],i] + term
                        
                    }
                }
            }

            CN <- c()

            for( f in 1:FN ){
                CN <- c(CN,colnames(OO[[f]]))
            }

            RES <- matrix(0,Cmax,length(CN))
            colnames(RES) <- CN

            RAW <- matrix(0,Cmax,length(CN))
            colnames(RAW) <- CN            
            
            for( f in 1:FN ){
                Cindx = match(colnames(OO[[f]]),colnames(RES))
                RES[,Cindx] = OO[[f]]
                RAW[,Cindx] = RW[[f]]
            }

            if( normROW ){
            
                for( i in 1:Cmax ){
                    Norm = sum(RES[i,])
                    if( Norm != 0 ){
                        RES[i,] = RES[i,]/Norm
                    }
                }
            }

        }        
    }

    return(list(RES=RES,RAW=RAW))
    
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

plotCor <- function(MAT=NULL, TARGETS=NULL, DIR="", subDIR="",
                    TIT="", shortTIT="", gradientTIT="Pearson\nCorrelation",
                    WIDTH=480, HEIGHT=480, valSIZE=3){

    gplot1   = NULL
    BORDERCOL= "grey50"#"black"#grey60
    LINECOL  = BORDERCOL
    LABELCOL = "black"
    NACOL    = "grey90"
    YMIN     = -4.0
    YMAX     = -4.0
    
    if( !is.null(MAT) ){

        CN = colnames(MAT)
        N  = length(CN)
        ii = seq_along(1:N)
        
        if( !is.null(TARGETS) ){

            Nt  = length(TARGETS)
            MIN = rep(NA,Nt)
            MAX = rep(NA,Nt)
            MID = rep(NA,Nt)

            LABELS = list()
            
            for( i in 1:Nt ){
                temp   = ii[grepl(TARGETS[i],CN)]
                MIN[i] = temp[1]
                MAX[i] = temp[length(temp)]
                MID[i] = MIN[i] + ((MAX[i]-MIN[i])/2)
                #cat(MIN[i], ", ", MAX[i], ", ", MID[i], "\n")
                LABELS[[i]] = textGrob(TARGETS[i], gp=gpar(fontsize=15, fontface="bold",col=LABELCOL))
            }
            
        }       

        Nlabs = length(LABELS)
        
        colnames(MAT) = ii
        rownames(MAT) = ii
                
        cormat    <- round(MAT,2)
   
    
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(cormat)#, na.rm=T)
    df = melted_cormat    
    colnames(df) <- c("V1","V2","value")
    
    # Create a ggheatmap
    gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
        labs(x="",y="",title=TIT)+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name=gradientTIT,#"Pearson\nCorrelation",
                             na.value=NACOL) +
        #theme_minimal()+ # minimal theme
        theme_grey(base_size=10)+
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
        coord_fixed()+        

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = BORDERCOL, fill=NA, size=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = BORDERCOL),
        axis.ticks = element_blank(),
        legend.justification = c("right","center"),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.position = "right",
        legend.direction = "vertical")+
    
        guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                     title.position = "top", title.hjust = 0.5))+
    
    {if( !is.na(MAX[1]))geom_vline(xintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[1]))geom_hline(yintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_vline(xintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_hline(yintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_vline(xintercept = MAX[3], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_hline(yintercept = MAX[3], color=LINECOL);}+ 
    {if( 1 <= Nlabs )annotation_custom(LABELS[[1]],xmin=MID[1],xmax=MID[1],ymin=YMIN,ymax=YMAX);}+
    {if( 2 <= Nlabs )annotation_custom(LABELS[[2]],xmin=MID[2],xmax=MID[2],ymin=YMIN,ymax=YMAX);}+
    {if( 3 <= Nlabs )annotation_custom(LABELS[[3]],xmin=MID[3],xmax=MID[3],ymin=YMIN,ymax=YMAX);}+
    {if( 4 <= Nlabs )annotation_custom(LABELS[[4]],xmin=MID[4],xmax=MID[4],ymin=YMIN,ymax=YMAX);}
         
        
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot1)
    dev.off()

    }
        
    return(gplot1)
}

plotMI <- function(MAT=NULL, TARGETS=NULL, DIR="", subDIR="",
                    TIT="", shortTIT="",
                    WIDTH=480, HEIGHT=480, valSIZE=3){

    gplot1   = NULL
    BORDERCOL= "grey50"#"black"#grey60
    LINECOL  = BORDERCOL
    LABELCOL = "black"
    NACOL    = "grey90"
    YMIN     = -4.0
    YMAX     = -4.0
    
    if( !is.null(MAT) ){

        CN = colnames(MAT)
        N  = length(CN)
        ii = seq_along(1:N)
        
        if( !is.null(TARGETS) ){

            Nt  = length(TARGETS)
            MIN = rep(NA,Nt)
            MAX = rep(NA,Nt)
            MID = rep(NA,Nt)

            LABELS = list()
            
            for( i in 1:Nt ){
                temp   = ii[grepl(TARGETS[i],CN)]
                MIN[i] = temp[1]
                MAX[i] = temp[length(temp)]
                MID[i] = MIN[i] + ((MAX[i]-MIN[i])/2)
                #cat(MIN[i], ", ", MAX[i], ", ", MID[i], "\n")
                LABELS[[i]] = textGrob(TARGETS[i], gp=gpar(fontsize=15, fontface="bold",col=LABELCOL))
            }
            
        }       

        Nlabs <- length(LABELS)
        
        colnames(MAT) = ii
        rownames(MAT) = ii
                
        cormat    <- round(MAT,2)
   
    
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(cormat)#, na.rm=T)
    df = melted_cormat    
    colnames(df) <- c("V1","V2","value")

    top = (max(as.numeric(as.vector(df$value))))
    bot = (min(as.numeric(as.vector(df$value))))
    mid = ((top-bot)/2) - abs(bot)
        
    # Create a ggheatmap
    gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
        labs(x="",y="",title=TIT)+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "cyan", high = "red", mid="black", 
                             midpoint = mid, limit = c(bot,top), space = "Lab", 
                             name="MI",
                             na.value=NACOL) +
        #theme_minimal()+ # minimal theme
        theme_grey(base_size=10)+
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
        coord_fixed()+        

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = BORDERCOL, fill=NA, size=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = BORDERCOL),
        axis.ticks = element_blank(),
        legend.justification = c("right","center"),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.position = "right",
        legend.direction = "vertical")+
    
        guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                     title.position = "top", title.hjust = 0.5))+
    
    {if( !is.na(MAX[1]))geom_vline(xintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[1]))geom_hline(yintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_vline(xintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_hline(yintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_vline(xintercept = MAX[3], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_hline(yintercept = MAX[3], color=LINECOL);}+ 
     {if( 1 <= Nlabs )annotation_custom(LABELS[[1]],xmin=MID[1],xmax=MID[1],ymin=YMIN,ymax=YMAX);}+
    {if( 2 <= Nlabs )annotation_custom(LABELS[[2]],xmin=MID[2],xmax=MID[2],ymin=YMIN,ymax=YMAX);}+
    {if( 3 <= Nlabs )annotation_custom(LABELS[[3]],xmin=MID[3],xmax=MID[3],ymin=YMIN,ymax=YMAX);}+
    {if( 4 <= Nlabs )annotation_custom(LABELS[[4]],xmin=MID[4],xmax=MID[4],ymin=YMIN,ymax=YMAX);}
         
        
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot1)
    dev.off()

    }
        
    return(gplot1)
}

plotII <- function(MAT=NULL, TARGETS=NULL, DIR="", subDIR="",
                    TIT="", shortTIT="",
                    WIDTH=480, HEIGHT=480, valSIZE=3){

    gplot1   = NULL
    BORDERCOL= "grey50"#"black"#grey60
    LINECOL  = BORDERCOL
    LABELCOL = "black"
    NACOL    = "grey90"
    YMIN     = -4.0
    YMAX     = -4.0
    
    if( !is.null(MAT) ){

        CN = colnames(MAT)
        N  = length(CN)
        ii = seq_along(1:N)
        
        if( !is.null(TARGETS) ){

            Nt  = length(TARGETS)
            MIN = rep(NA,Nt)
            MAX = rep(NA,Nt)
            MID = rep(NA,Nt)

            LABELS = list()
            
            for( i in 1:Nt ){
                temp   = ii[grepl(TARGETS[i],CN)]
                MIN[i] = temp[1]
                MAX[i] = temp[length(temp)]
                MID[i] = MIN[i] + ((MAX[i]-MIN[i])/2)
                #cat(MIN[i], ", ", MAX[i], ", ", MID[i], "\n")
                LABELS[[i]] = textGrob(TARGETS[i], gp=gpar(fontsize=15, fontface="bold",col=LABELCOL))
            }
            
        }       

        Nlabs <- length(LABELS)
        
        colnames(MAT) = ii
        rownames(MAT) = ii
                
        cormat    <- round(MAT,2)
   
    
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(cormat)#, na.rm=T)
    df = melted_cormat    
    colnames(df) <- c("V1","V2","value")

    top = (max(as.numeric(as.vector(df$value))))
    bot = (min(as.numeric(as.vector(df$value))))
    mid = ((top-bot)/2) - abs(bot)
        
    # Create a ggheatmap
    gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
        labs(x="",y="",title=TIT)+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "cyan", high = "red",  
                             limit = c(bot,top), space = "Lab", 
                             name="II",
                             na.value=NACOL) +
        #theme_minimal()+ # minimal theme
        theme_grey(base_size=10)+
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
        coord_fixed()+        

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = BORDERCOL, fill=NA, size=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = BORDERCOL),
        axis.ticks = element_blank(),
        legend.justification = c("right","center"),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.position = "right",
        legend.direction = "vertical")+
    
        guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                     title.position = "top", title.hjust = 0.5))+
    
    {if( !is.na(MAX[1]))geom_vline(xintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[1]))geom_hline(yintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_vline(xintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_hline(yintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_vline(xintercept = MAX[3], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_hline(yintercept = MAX[3], color=LINECOL);}+ 
     {if( 1 <= Nlabs )annotation_custom(LABELS[[1]],xmin=MID[1],xmax=MID[1],ymin=YMIN,ymax=YMAX);}+
    {if( 2 <= Nlabs )annotation_custom(LABELS[[2]],xmin=MID[2],xmax=MID[2],ymin=YMIN,ymax=YMAX);}+
    {if( 3 <= Nlabs )annotation_custom(LABELS[[3]],xmin=MID[3],xmax=MID[3],ymin=YMIN,ymax=YMAX);}+
    {if( 4 <= Nlabs )annotation_custom(LABELS[[4]],xmin=MID[4],xmax=MID[4],ymin=YMIN,ymax=YMAX);}
    
       
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot1)
    dev.off()

    }
        
    return(gplot1)
}


plotNMI <- function(MAT=NULL, TARGETS=NULL, DIR="", subDIR="",
                    TIT="", shortTIT="", gradientTIT="NMI",
                    WIDTH=480, HEIGHT=480, valSIZE=3){

    gplot1   = NULL
    BORDERCOL= "grey50"#"black"#grey60
    LINECOL  = BORDERCOL
    LABELCOL = "black"
    NACOL    = "grey90"
    YMIN     = -4.0
    YMAX     = -4.0
    
    if( !is.null(MAT) ){

        CN = colnames(MAT)
        N  = length(CN)
        ii = seq_along(1:N)
        
        if( !is.null(TARGETS) ){

            Nt  = length(TARGETS)
            MIN = rep(NA,Nt)
            MAX = rep(NA,Nt)
            MID = rep(NA,Nt)

            LABELS = list()
            
            for( i in 1:Nt ){
                temp   = ii[grepl(TARGETS[i],CN)]
                MIN[i] = temp[1]
                MAX[i] = temp[length(temp)]
                MID[i] = MIN[i] + ((MAX[i]-MIN[i])/2)
                #cat(MIN[i], ", ", MAX[i], ", ", MID[i], "\n")
                LABELS[[i]] = textGrob(TARGETS[i], gp=gpar(fontsize=15, fontface="bold",col=LABELCOL))
            }
            
        }       

        Nlabs = length(LABELS)
        
        colnames(MAT) = ii
        rownames(MAT) = ii
                
        cormat    <- round(MAT,2)
   
    
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(cormat)#, na.rm=T)
    df = melted_cormat    
    colnames(df) <- c("V1","V2","value")
    
    # Create a ggheatmap
    gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
        labs(x="",y="",title=TIT)+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "cyan", high = "red", mid="black", 
                             midpoint = 0.5, limit = c(0,1), space = "Lab", 
                             name=gradientTIT,
                             na.value=NACOL) +
        #theme_minimal()+ # minimal theme
        theme_grey(base_size=10)+
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
        coord_fixed()+        

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = BORDERCOL, fill=NA, size=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = BORDERCOL),
        axis.ticks = element_blank(),
        legend.justification = c("right","center"),
        legend.margin=margin(grid::unit(0,"cm")),
        legend.position = "right",
        legend.direction = "vertical")+
    
        guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                                     title.position = "top", title.hjust = 0.5))+
    
    {if( !is.na(MAX[1]))geom_vline(xintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[1]))geom_hline(yintercept = MAX[1], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_vline(xintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[2]))geom_hline(yintercept = MAX[2], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_vline(xintercept = MAX[3], color=LINECOL);}+
    {if( !is.na(MAX[3]))geom_hline(yintercept = MAX[3], color=LINECOL);}+ 
      {if( 1 <= Nlabs )annotation_custom(LABELS[[1]],xmin=MID[1],xmax=MID[1],ymin=YMIN,ymax=YMAX);}+
    {if( 2 <= Nlabs )annotation_custom(LABELS[[2]],xmin=MID[2],xmax=MID[2],ymin=YMIN,ymax=YMAX);}+
    {if( 3 <= Nlabs )annotation_custom(LABELS[[3]],xmin=MID[3],xmax=MID[3],ymin=YMIN,ymax=YMAX);}+
    {if( 4 <= Nlabs )annotation_custom(LABELS[[4]],xmin=MID[4],xmax=MID[4],ymin=YMIN,ymax=YMAX);}
        
        
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot1)
    dev.off()

    }
        
    return(gplot1)
}


binomialModel <- function( data, alpha0=NULL, beta0=NULL, plot=TRUE ){
##Ref:  http://dm13450.github.io/2019/02/22/Nonparametric-Prior.html    

    gplot=NULL
    
    if( is.null(alpha0) ){ alpha0 = 0.01 }
    if( is.null(beta0)  ){ beta0  = 0.01 }

    alphaPosterior <- alpha0 + sum(data[,1]) 
    betaPosterior  <- beta0  + sum(data[,2]) - sum(data[,1])

    thetaDraws <- rbeta(1000, alphaPosterior, betaPosterior)

    if( plot ){
    
        xGrid <- seq(0.01, 1-0.01, 0.01)
    
        bayesFrame <- data.frame(Theta=xGrid,
                                 Prior = dbeta(xGrid, alpha0, beta0),
                                 Likelihood = vapply(xGrid, function(x) sum(dbinom(data[,1], data[,2], x, log=T)), numeric(1)),
                                 Posterior = vapply(xGrid, function(x) dbeta(x, alphaPosterior, betaPosterior), numeric(1)))

        bayesFrame %>% gather(Part, Value, -Theta) -> bayesFrameTidy
        bayesFrameTidy$Part <- factor(bayesFrameTidy$Part, levels=c("Prior", "Likelihood", "Posterior"))

        gplot = ggplot(bayesFrameTidy, aes(x=Theta, y=Value, colour=Part)) + 
            geom_line() + 
            facet_wrap(~Part, scales="free") + 
            theme(legend.position = "bottom")

    }

    return(list(THETA=thetaDraws,GPLOT=gplot))
    
}

plotJcmi <- function(DF, DIR="", subDIR="", TIT="", shortTIT="", WIDTH=480, HEIGHT=480, STEPS=0.05, INTmax=30, runSCALED=FALSE, scaledMIN=-0.5, scaledMAX=0.5 ){

    DF <- as.data.frame(DF)

    CLASSES <- levels(DF[,2])

    PLOTS <- list()
    
    for( i in 1:length(CLASSES) ){

        DF2 = DF[grepl(CLASSES[i], DF[,2]),]
        
        MIN <- min(as.numeric(as.vector(DF2[,1])))
        MAX <- max(as.numeric(as.vector(DF2[,1])))
        INT <- ceiling((MAX-MIN)/STEPS)

        #print(INT)
        
        if( INT > INTmax ){ INT = INTmax }
        
        BW  <- round((MAX-MIN)/INT,3)    

        gplot <- ggplot(DF2, aes(as.numeric(as.vector(jcmi))))+#,fill=class))+
            geom_histogram(binwidth=BW,fill="gray80",color="black")+    
        scale_x_continuous(breaks = seq(MIN, MAX, by = BW))+
        stat_bin(binwidth = BW, aes(label=..count..), vjust=-0.5, geom = "text")+
            geom_vline(xintercept = 0, color = "red", size=1)+
    labs(x=shortTIT,y="count",title=CLASSES[i])+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.0)),
            axis.title.y=element_text(face="bold",size=rel(1.0)),
            axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
            legend.title=element_text(face="bold",size=rel(1.0)),
            legend.text=element_text(face="bold",size=rel(1.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              #panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))
        
        png(sprintf("%s/%s/%s_%s.png",DIR,subDIR,shortTIT,CLASSES[i]),width=WIDTH,height=HEIGHT,units="px");
        print(gplot)
        dev.off()

        PLOTS[[i]] <- gplot
        names(PLOTS)[i] <- CLASSES[i]
        
    }

    MIN <- min(as.numeric(as.vector(DF[,1])))
    MAX <- max(as.numeric(as.vector(DF[,1])))
    
    gplot <- ggplot(DF) +
        geom_density(aes(x = as.numeric(as.vector(jcmi)), fill = class), alpha = 0.2)+
        geom_vline(xintercept = 0, color = "red", size=1)+
        labs(x=shortTIT,y="count",title="ALL Classes")+
        xlim(c(MIN,MAX))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.0)),
            axis.title.y=element_text(face="bold",size=rel(1.0)),
            axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
            legend.title=element_text(face="bold",size=rel(1.0)),
            legend.text=element_text(face="bold",size=rel(1.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              #panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))
        
    png(sprintf("%s/%s/%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot)
    dev.off()

    i=(i+1)
    
    PLOTS[[i]]      <- gplot
    names(PLOTS)[i] <- "ALL"

    if( runSCALED ){

        MIN = scaledMIN
        MAX = scaledMAX
        
        gplot <- ggplot(DF) +
            geom_density(aes(x = as.numeric(as.vector(jcmi)), fill = class), alpha = 0.2)+
            geom_vline(xintercept = 0, color = "red", size=1)+
            labs(x=shortTIT,y="count",title="ALL Classes")+
           xlim(c(MIN,MAX))+
            theme(
                title=element_text(face="bold",size=rel(1.5)),
                axis.title.x=element_text(face="bold",size=rel(1.0)),
                axis.title.y=element_text(face="bold",size=rel(1.0)),
                axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.0)),
                legend.title=element_text(face="bold",size=rel(1.0)),
                legend.text=element_text(face="bold",size=rel(1.0)),
                legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              #panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_rect(linetype="solid",fill=NA))
        
        png(sprintf("%s/%s/%s_scaled.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
        print(gplot)
        dev.off()

        i=(i+1)
    
        PLOTS[[i]]      <- gplot
        names(PLOTS)[i] <- "ALLscaled"
        
        
    }
    
    return(PLOTS)
    
}

binomialDPModel <- function( data, alpha0=NULL, beta0=NULL, ITS=20, FIT=TRUE ){
##Ref:  http://dm13450.github.io/2019/02/22/Nonparametric-Prior.html    

    saveList <- NULL
    
    if(FIT){

        xGrid <- seq(0,1,by=0.01)
        
        if( is.null(alpha0) ){ alpha0 = 0.01 }
        if( is.null(beta0)  ){ beta0  = 0.01 }

        N = length(data[,1])
        
        alphaPosterior <- alpha0 + sum(data[,1]) 
        betaPosterior  <- beta0  + sum(data[,2]) - sum(data[,1])
    
        thetaDirichlet <- rbeta(N, alphaPosterior, betaPosterior)

        dp <- DirichletProcessBeta( thetaDirichlet, 
                                   1,
                                   mhStep = c(0.002, 0.005),
                                   #(alpha,beta) shape parameters for beta distribution
                                   #https://en.wikipedia.org/wiki/Beta_distribution
                                   alphaPrior = c(2,2))#c(2,0.5) )
        
        dp$numberClusters    <- N
        dpclusterLabels      <- seq_len(N)
        dp$clusterParameters <- PriorDraw(dp$mixingDistribution, N)
        dp$pointsPerCluster  <- rep_len(1, N)

        dp                   <- Fit(dp,1)

        postFuncEval <- matrix(ncol=ITS, nrow=length(xGrid))
        muPostVals   <- matrix(ncol=ITS, nrow=N)
        nuPostVals   <- matrix(ncol=ITS, nrow=N)

        pb <- txtProgressBar(max=ITS, width=50, char="-", style=3)

        Error   = 1
        Problem = PosteriorClusters(dp)$weights
        Theta   = thetaDirichlet
        
        tryCatch(

            expr = {
                for( i in seq_len(ITS) ){
                    postClusters     <- PosteriorClusters(dp)
                    postFuncEval[,i] <- PosteriorFunction(dp)(xGrid)

                    if( length(postClusters$weights) > 0 &&
                        sum(is.na(postClusters$weights)) == 0 ){
                    
                        wk <- sample.int(length(postClusters$weights),
                                         N,
                                         replace=T,
                                         prob=postClusters$weights)
                        
                        muPost           <- postClusters$params[[1]][,,wk]
                        nuPost           <- postClusters$params[[2]][,,wk]
                        
                        aPost            <- muPost * nuPost
                        bPost            <- (1-muPost) * nuPost
                        
                        muPostVals[,i]   <- muPost
                        nuPostVals[,i]   <- nuPost

                        newTheta         <- rbeta(N, aPost+data[,1], bPost+data[,2] - data[,1])
                        dp               <- ChangeObservations(dp, newTheta)
                        dp               <- Fit(dp, 100, updatePrior=T, progressBar=F)
                    
                        Error            <- 0
                        Problem          <- postClusters$weights
                        Theta            <- newTheta

                        #generate an error
                        #log("text")
                        
                    } 
                        
                    setTxtProgressBar(pb, i)
                }#for
                
            },

            error = function(e){
                message('Caught an error!')
                print(e)
                Error = 1
                Problem = postClusters$weights
            }
           #,
           #
           # warning = function(w){
           #     message('Caught an warning!')
           #     #print(w)
           # }                        
       )        
        
        saveList              <- list()
        saveList$muPostVals   <- muPostVals
        saveList$nuPostVals   <- nuPostVals
        saveList$postFuncEval <- postFuncEval
        saveList$dp           <- dp
        saveList$Error        <- Error
        saveList$Problem      <- Problem
        saveList$Theta        <- Theta
    }    

    return(saveList)
    
}


variablePlotBiomModel <- function(DF, DIR, subDIR, WIDTH, HEIGTH ){

    gplot = ggplot(DF, aes(x=as.numeric(as.vector(mean))) ) +
        stat_bin(breaks=c(-0.004, seq(0.001,1.0, by=0.005)),color="#E69F00")+
        labs(x="mean(theta)",y="count")+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.5)),
            axis.title.y=element_text(face="bold",size=rel(1.5)),
            legend.title=element_text(face="bold",size=rel(2.0)),
            legend.text=element_text(face="bold",size=rel(2.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))

    png(sprintf("%s/%s/%s.png",DIR,subDIR,"Mn_theta"),width=WIDTH,height=HEIGHT,units="px");
    print(gplot)
    dev.off()


    gplot = ggplot(DF, aes(x=as.numeric(as.vector(mean)),color=set, fill=set) ) +
        stat_bin(breaks=c(-0.004, seq(0.001,1.0, by=0.005)))+
        labs(x="mean(theta)",y="count")+
        scale_color_manual(values=c("#E69F00", "#999999", "#56B4E9", "#CC79A7"))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(1.5)),
            axis.title.y=element_text(face="bold",size=rel(1.5)),
            legend.title=element_text(face="bold",size=rel(2.0)),
            legend.text=element_text(face="bold",size=rel(2.0)),
            legend.key=element_blank(),
            legend.background=element_blank() )+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))

    png(sprintf("%s/%s/%s.png",DIR,subDIR,"Mn_theta_bySet"),width=WIDTH,height=HEIGHT,units="px");
    print(gplot)
    dev.off()
    
}
