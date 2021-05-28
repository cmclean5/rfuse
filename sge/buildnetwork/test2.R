## run setParameters.R
source("setParameters.R")

## run setUp.R
source("setUp.R")

#----------------------------------
# compile c++ code (fuse.cpp)
#----------------------------------

#Sys.setenv("PKG_CXXFLAGS" = cxxflag.opts)
#Sys.setenv("PKG_LIBS"     = lib.opts)
#sourceCpp(sprintf("%s/%s",PATHS[4],files[6]))

library(rfuse)
library(Matrix) #for forceSymmetric

#create tmp dir to hold matrices
td <- tempdir()
print(td)

NR=4
NR2=NR*NR

tmp1 = "tmp"
st1 = sprintf("%s/%s.bk",td,tmp1)
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }
st1 = sprintf("%s/%s.rds",td,tmp1)
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }

tmp = FBM(NR,NR,backingfile=file.path(td,tmp1),"double")

tmp2 = "tmp2"
st1 = sprintf("%s/%s.bk",td,tmp2)
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }
st1 = sprintf("%s/%s.rds",td,tmp2)
if( file.exists(st1) && file.info(st1)$size!=0 ){ file.remove(st1); }

tmp2 = FBM(NR,NR,backingfile=file.path(td,tmp2),"double")


x=rnorm(NR2)
X=matrix(x,NR,NR)
X=forceSymmetric(X)

X2=X
X2[sample(NR2,10)] = 0

#X  = matrix(1,NR,NR)
#X2 = X
#X2[sample(NR2,10)] = 0

print(X)
print(X2)

for( i in 1:NR ){
    for(j in 1:NR ){
        tmp[i,j]  = X[i,j]
        tmp2[i,j] = X2[i,j]
    }
}

createModel()

test(tmp,tmp2,c(1),c(2))

deleteModel()

z = t(X2) %*% X %*% X2

cat("is t(X2) * X * X2 symmetric: ")
if( Matrix::isSymmetric(z) == TRUE ){
    cat("YES")
} else {
    cat("NO")
}
cat("\n")

#test4(tmp,tmp2,c(1),c(2))


#unlink(sprintf("%s/%s",file.path(td),paste0("tmp" ,c(".bk",".rds"))))
