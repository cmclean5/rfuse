//include external header's gloabl parameters, and c/c++ functions
#include "fuse.h"

fuse::fuse(){

  //variables
  EPSILON    = std::numeric_limits<double>::min();
  logEPSILON = log(std::numeric_limits<double>::min());

  OMP_start  = 0;
  OMP_end    = 0;

  tSr1csr    = nullptr;
  Sr1csc     = nullptr;
  tSr2csr    = nullptr;
  Sr2csc     = nullptr;
  tSr3csr    = nullptr;
  Sr3csc     = nullptr;
  tSr4csr    = nullptr;
  Sr4csc     = nullptr;
  
  PrDiag1    = nullptr;
  PrDiag2    = nullptr;
  PrDiag3    = nullptr;
  PrDiag4    = nullptr;

  PrBlck1    = nullptr;
  PrBlck2    = nullptr;
  PrBlck3    = nullptr;
  PrBlck4    = nullptr;

  freedCSFSpace=false;
  freedPrDiagSpace=false;
  freedPrBlckSpace=false;
  
}

fuse::~fuse(){

  freeSpace();

}

void fuse::setOpenMP( int ncores, bool print ){

  std::unordered_map<unsigned,std::string> map{
    {200505,"2.5"},{200805,"3.0"},{201107,"3.1"},{201307,"4.0"},{201511,"4.5"},{201811,"5.0"},{202011,"5.1"}};
  if( print ){ cout << "--------------" << endl; }
  if( print ){ cout << "OpenMP: V" << map.at(_OPENMP) << ". "; }

  omp_set_dynamic(false);
  
  //int ncores      = nCORES[0];
  int max_threads = 0;
  int num_procs   = 0;
  max_threads     = omp_get_max_threads();
  num_procs       = omp_get_num_procs();  
  if( (ncores < 0) || (ncores > max_threads) ){
    omp_set_num_threads(max_threads);
    if( print){ cout << "cores set too: " << max_threads << endl; }
  } else {
    if( ncores == 0 ) ncores = 1;
    omp_set_num_threads(ncores);
    if( print ){ cout << "cores set too: " << ncores << endl; }
    if( print ){ cout << "--------------" << endl; }
  }
  //---  

 
}

void fuse::getOpenMPStart( double &omp_start ){  omp_start = omp_get_wtime(); }

void fuse::getOpenMPEnd  ( double &omp_end   ){  omp_end   = omp_get_wtime(); }


void fuse::freeCSFSpace(){

  //need to redo - hiding the fact that i'm double deleting. 
  if(!freedCSFSpace){
  
  if(tSr1csr!=nullptr){ delete tSr1csr; }//tSr1csr=nullptr; }
  if( Sr1csc!=nullptr){ delete  Sr1csc; }//Sr1csc=nullptr;  }
 
  if(tSr2csr!=nullptr){ delete tSr2csr; }//tSr2csr=nullptr; }
  if( Sr2csc!=nullptr){ delete  Sr2csc; }//Sr2csc=nullptr;  }
   
  if(tSr3csr!=nullptr){ delete tSr3csr; }//tSr3csr=nullptr; }
  if( Sr3csc!=nullptr){ delete  Sr3csc; }//Sr3csc=nullptr;  }
 
  if(tSr4csr!=nullptr){ delete tSr4csr; }//tSr4csr=nullptr; }
  if( Sr4csc!=nullptr){ delete  Sr4csc; }//Sr4csc=nullptr;  }

  freedCSFSpace=true;
  
  }
  
}


void fuse::freePrDiagSpace(){

  if(!freedPrDiagSpace){
  
  if(PrDiag1!=nullptr){free(PrDiag1); }//PrDiag1=nullptr; }
  if(PrDiag2!=nullptr){free(PrDiag2); }//PrDiag2=nullptr; }
  if(PrDiag3!=nullptr){free(PrDiag3); }//PrDiag3=nullptr; }
  if(PrDiag4!=nullptr){free(PrDiag4); }//PrDiag4=nullptr; }

  freedPrDiagSpace=true;
  
  }
  
}

void fuse::freePrBlckSpace(){

  if(!freedPrBlckSpace){
  
  if(PrBlck1!=nullptr){_mm_free(PrBlck1); }//PrBlck1=nullptr; }
  if(PrBlck2!=nullptr){_mm_free(PrBlck2); }//PrBlck2=nullptr; }
  if(PrBlck3!=nullptr){_mm_free(PrBlck3); }//PrBlck3=nullptr; }
  if(PrBlck4!=nullptr){_mm_free(PrBlck4); }//PrBlck4=nullptr; }

  freedPrBlckSpace=true;
  
  }
  
}

void fuse::freeSrSpace(){
 
  
  if(mapSr1.capacity()>0){std::vector<tupleSr>().swap(mapSr1);}
  if(mapSr2.capacity()>0){std::vector<tupleSr>().swap(mapSr2);}
  if(mapSr3.capacity()>0){std::vector<tupleSr>().swap(mapSr3);}
  if(mapSr4.capacity()>0){std::vector<tupleSr>().swap(mapSr4);}
  
  
}

void fuse::freeSrCOOSpace(){

  //free Sr memory vectors
  if(COOSr1.capacity()>0){std::vector<tupleCOO>().swap(COOSr1);}
  if(COOSr2.capacity()>0){std::vector<tupleCOO>().swap(COOSr2);}
  if(COOSr3.capacity()>0){std::vector<tupleCOO>().swap(COOSr3);}
  if(COOSr4.capacity()>0){std::vector<tupleCOO>().swap(COOSr4);}
  
  
}


void fuse::freeMapSpace(){

  //free map space
  if(map2Arr.capacity()>0){std::vector<tripleInt>().swap(map2Arr);}
  
}

void fuse::setMapSpace( int K ){

  map2Arr.reserve(K);
  
}

void fuse::freeSpace(){

  freeCSFSpace();
  freePrDiagSpace();
  freePrBlckSpace();
  freeSrSpace();
  freeSrCOOSpace();
  freeMapSpace();
   

}


void fuse::copyDatasetRtoC( double *Dnew, BMAcc_RW<double> &Dold ){

  //int i,N,K;
  uint_fast32_t i,N,K;
  
  N = Dold.nrow();
  K = N*N;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i) shared(Dold,Dnew)
  for(i=0; i<K; i++){
    Dnew[i] = 0;
    Dnew[i] = Dold[i];
  }
  
  
}


void fuse::copyDatasetRtoC( BMAcc_RW<double> &Dold, double *Pnew, const uint_fast32_t n, const vector<tripleInt> &map, double &sum, bool init ){

  uint_fast32_t i,N,K;
  
  N = n;
  K = N*N;  
  
  double *Dnew = (double*)malloc((N*N)*sizeof(double));
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i) shared(Dold,Dnew)
  for(i=0; i<K; i++){
    Dnew[i] = 0;
    Dnew[i] = Dold[i];
  }

  if( init == true ){

    cout << " (initialised) "; 
    
    //row normalisation
    rowNorm( Dnew, N );

    //symmetrise
    symmetrise( Dnew, N, sum );

  }
  
  map2Array( Dnew, Pnew, N, map );

  //free space
  if(Dnew){free(Dnew);}
  //if(Dnew){_mm_free(Dnew);}
  
}

void fuse::csr_dense_tcrossprod( const CSF *Xcsr, const arma::Mat<double> &tY,
				 arma::Mat<double> &res, const int nTHREADS ){
 //https://github.com/rexyai/rsparse/blob/a887082b52f4fbe5d8b22a0b19c8d2d07e630d6e/src/matrix_product.cpp
  // implements X * t(Y) multiplication where X = sparse CSR matrix and Y = dense matrix
  // it doesn't implement just X * Y because it will be slow due to the fact that
  // our main pattern for access elements in Y is by row which is cache inefficient and
  // problematic to vectorize because matrices are stored by column so essentially it
  // implements R's `tcrossprod()`, i.e. t(tcrossprod(Y,X))
  // The equivalent using Armadillo is: X * Y

  arma::uword ii, NR;

  NR = Xcsr->nrows;
  
  res.zeros(Xcsr->nrows,tY.n_rows);
  
  //#pragma omp parallel for schedule(static) default(none) firstprivate(NR) private(ii) shared(Xcsr,tY,res)
#pragma omp parallel for num_threads(nTHREADS) schedule(dynamic, CHUNK) default(none) firstprivate(NR) private(ii) shared(Xcsr,tY,res)
  for(ii=0; ii<NR; ii++) {
    const arma::uword p1 = Xcsr->row[ii];
    const arma::uword p2 = Xcsr->row[ii+1];
    const arma::uvec idx = arma::uvec(&Xcsr->col[p1], p2 - p1, false, true);
    const arma::colvec csr_row = arma::colvec(&Xcsr->val[p1], p2 - p1, false, false);
    res.row(ii) = (tY.cols(idx) * csr_row).t();
  }  
  
}

void fuse::dense_csc_prod(const arma::Mat<double>& X, const CSF *Ycsc,
			  arma::Mat<double> &res, const int nTHREADS ){

  arma::uword ii, NC;

  NC = Ycsc->ncols;

  res.zeros(X.n_rows,Ycsc->ncols); 
  
  //#pragma omp parallel for schedule(static) default(none) firstprivate(NC) private(ii) shared(X,Ycsc,res)
  #pragma omp parallel for num_threads(nTHREADS) schedule(dynamic, CHUNK) default(none) firstprivate(NC) private(ii) shared(X,Ycsc,res)
  for (ii=0; ii<NC; ii++) {
    const arma::uword p1 = Ycsc->col[ii];
    const arma::uword p2 = Ycsc->col[ii+1];
    const arma::uvec idx = arma::uvec(&Ycsc->row[p1], p2 - p1, false, true);
    const arma::colvec csc_col = arma::colvec(&Ycsc->val[p1], p2 - p1, false, false);
    res.col(ii) = X.cols(idx) * csc_col;
  }
  

}

void fuse::mappingCOO2CSF( vector<tupleCOO> &COOmap, CSF *CSFmap, const uint_fast32_t n, bool colMajor ){

  arma::uword i,k,COOsize,tally,tmp,N,NP1;

  N   = static_cast<arma::uword>(n);
  NP1 = static_cast<arma::uword>(N+1); 
  COOsize = static_cast<arma::uword>(COOmap.size()); 
   
  
  if( colMajor ){ //Compressed Column Storage format   

    //order column indices
    sort(COOmap.begin(),COOmap.end(),sortColCOO());
      
    //initialise CSF struct
    CSFmap->init_csc( COOsize, N );   

    #pragma omp parallel for schedule(static) default(none) firstprivate(COOsize) private(k) shared(COOmap,CSFmap)
    for(k=0; k<COOsize; k++){
      CSFmap->val[k] = std::get<2>(COOmap[k]);
      CSFmap->row[k] = std::get<0>(COOmap[k]);    
    }

   
    #pragma omp parallel for schedule(static) default(none) private(i) firstprivate(NP1) shared(CSFmap)
    for(i=0;i<NP1; i++){
      CSFmap->col[i]=0; 
    }
    
    #pragma omp parallel for schedule(static) default(none) private(i,k) firstprivate(COOsize) shared(COOmap,CSFmap)
    for(k=0;k<COOsize;k++){
      i = std::get<1>(COOmap[k]);
      CSFmap->col[i]++;
    }
    
    tmp    = 0;
    tally  = 0;

    ////#pragma omp parallel for schedule(static) default(none) private(i,tmp) firstprivate(NP1) shared(CSFmap) reduction(+:tally)
    for(i=0; i<NP1; i++){
      if( i == 0 ){
	tmp = 0;
	tally = CSFmap->col[i];
	CSFmap->col[i] = 0;
      } else {
	tmp = CSFmap->col[i];
	CSFmap->col[i] = tally;
	tally += tmp;
      }
    }

    
  } else { //Compressed Row Storage format

    //order row indices
    sort(COOmap.begin(),COOmap.end(),sortRowCOO());
    
    //initialise CSF struct
    CSFmap->init_csr( COOsize, N );    

    #pragma omp parallel for schedule(static) default(none) firstprivate(COOsize) private(k) shared(COOmap,CSFmap)
    for(k=0; k<COOsize; k++){
      CSFmap->val[k] = std::get<2>(COOmap[k]);
      CSFmap->col[k] = std::get<1>(COOmap[k]);    
    }
    

    #pragma omp parallel for schedule(static) default(none) private(i) firstprivate(NP1) shared(CSFmap)
    for(i=0;i<NP1; i++){
      CSFmap->row[i]=0; 
    }
    
    #pragma omp parallel for schedule(static) default(none) private(i,k) firstprivate(COOsize) shared(COOmap,CSFmap)
    for(k=0;k<COOsize;k++){
      i = std::get<0>(COOmap[k]);
      CSFmap->row[i]++;
    }
    
    tmp    = 0;
    tally  = 0;
    
    
    /////#pragma omp parallel for schedule(static) default(none) private(i,tmp) firstprivate(NP1) shared(CSFmap) reduction(storeTally:myTally)
    for(i=0; i<NP1; i++){
      if( i == 0 ){
      tmp = 0;
      tally = CSFmap->row[i];
      CSFmap->row[i] = 0;
      } else {
	tmp = CSFmap->row[i];
	CSFmap->row[i] = tally;
	tally += tmp;
      }
    }

  }//ifelse
 

}


void fuse::mapSrmatrix(BMAcc_RW<double> &Smat, const uint_fast32_t n, vector<tupleCOO> &map, const vector<tripleInt> &mapP, bool colMajor){
  
 uint_fast32_t i,j,k,kk,N,KK,KP,K,rind,cind;

 double Val, rVal;
  
 int total = 0;
  
  N  = n;
  K  = N*N;
  KK = N*(N+1)/2;
  KP = mapP.size();

   
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,Val) shared(Smat) reduction(+:total)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){

      rind = (i*N)+j;
      Val  = Smat[rind];

      if( Val != 0 ){ total++; }
    
    }
  }
  
  cout << "" << endl;
  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;
  
  map.reserve(total);

  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;

  // user-defined reduction for pushing back in vector map
  //https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
#pragma omp declare reduction (merge : std::vector<tupleCOO> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  if( colMajor ){ //column-major storage

    #pragma omp parallel for collapse(2) default(none) firstprivate(N) private(i,j,cind,Val) shared(Smat) reduction(merge: map)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	
	  cind = j*N+i; //column-major indexing, the transposeof row-amjor
	    
	  Val  = Smat[cind];	
	  if( Val != 0 ){ map.push_back( tupleCOO(i, j, Val) ); }
	 
      }
    }
    
  
  } else { //row-major storage

#pragma omp parallel for collapse(2) default(none) firstprivate(N) private(i,j,rind,Val) shared(Smat) reduction(merge: map)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){

	  rind = i*N+j; //because 'i' is our outer-loop, and 'j' our fast loop this is row-major indexing

	  Val = Smat[rind];	
	  if( Val != 0 ){ map.push_back( tupleCOO(i, j, Val) ); }
	
      }	  
    }
  }//ifelse 
    
  
  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;


}

void fuse::mapSrmatrix( BMAcc_RW<double> &Srmat, const uint_fast32_t N,
		  CSF *tCSFcsr, CSF *CSFcsc, const vector<tripleInt> &map ){

  vector<tupleCOO> Srmap, tSrmap;
  
  mapSrmatrix( Srmat, N, Srmap, map, true );//read Srmat in column-major format into COO format

  tSrVec( Srmap, tSrmap, N );//transpose Srmap

  mappingCOO2CSF( tSrmap, tCSFcsr, N, false );//store tSr in CSR format
   
  mappingCOO2CSF(  Srmap,  CSFcsc, N, true  );//store Sr in CSC format

  //free mapping vector
  std::vector<tupleCOO>().swap(Srmap);
  std::vector<tupleCOO>().swap(tSrmap);
  
}

void fuse::mapSrmatrix(BMAcc_RW<double> &Smat, const uint_fast32_t n, vector<tripleSr>&map, const vector<tripleInt> &mapP){
  
 uint_fast32_t i,j,k,kk,N,KK,KP,K,ind,rind;

  double Vind, Vrind;
  
  int total = 0;
  
  N  = n;
  K  = N*N;
  KK = N*(N+1)/2;
  KP = mapP.size();

   
#pragma omp parallel for schedule(static) default(none) firstprivate(N,KP) private(i,j,k,kk,ind,Vind) shared(Smat,mapP) reduction(+:total)
  for( kk=0; kk<KP; kk++ ){
    i = std::get<0>(mapP[kk]);
    j = std::get<1>(mapP[kk]);
    k = std::get<2>(mapP[kk]);

    ind  = (i*N)+j;    
    Vind = Smat[ind];     

    if( Vind != 0 ){ total++; }	
    
  }
  

  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;
  
  map.reserve(total);

  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;

  // user-defined reduction for pushing back in vector map
  //https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
#pragma omp declare reduction (merge : std::vector<tripleSr> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  #pragma omp parallel for schedule(static) default(none) firstprivate(N,KP) private(i,j,k,kk,ind,rind,Vind,Vrind) shared(Smat,mapP) reduction(merge: map)
  for( kk=0; kk<KP; kk++ ){
    i = std::get<0>(mapP[kk]);
    j = std::get<1>(mapP[kk]);
    k = std::get<2>(mapP[kk]);

    ind  = (i*N)+j;
    rind = (j*N)+i;

     Vind  = Smat[ind];
     Vrind = Smat[rind];

     if( Vind != 0 && Vrind != 0 ){
	
       map.push_back( tripleSr(k,Vind,Vrind) );
	    
     }
    
  }
 
  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;
    
 
}

void fuse::mapSrmatrix(BMAcc_RW<double> &Smat, const uint_fast32_t n, vector<tupleSr>&map, const vector<tripleInt> &mapP){

  uint_fast32_t i,j,k,kk,N,K,rind,cind;

  double total, Vrind, Vcind;
  
  N  = n;
  K  = mapP.size();  
  
  #pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(i,j,k,kk,rind,Vrind) shared(Smat,mapP) reduction(+:total)
  for( kk=0; kk<K; kk++ ){
    i = std::get<0>(mapP[kk]);
    j = std::get<1>(mapP[kk]);
    k = std::get<2>(mapP[kk]);

    rind  = (i*N)+j;    
    Vrind = Smat[rind];     

    if( Vrind != 0 ){ total++; }	
    
  }
  
  
  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;
  
  map.reserve(total);

  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;

  // user-defined reduction for pushing back in vector map
  //https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
#pragma omp declare reduction (merge : std::vector<tupleSr> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
   
  #pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(i,j,k,kk,rind,cind,Vrind,Vcind) shared(Smat,mapP) reduction(merge: map)
  for( kk=0; kk<K; kk++ ){
    i = std::get<0>(mapP[kk]);
    j = std::get<1>(mapP[kk]);
    k = std::get<2>(mapP[kk]);

    rind = (i*N)+j;
    cind = (j*N)+i;

    Vrind = Smat[rind];
    Vcind = Smat[cind];

    if( Vrind != 0 || Vcind != 0 ){    
      map.push_back( tupleSr(i,j,k,Vrind,Vcind) );
    }
    
  }     
 
  cout << "map capacity: " << map.capacity() << " map.size() " << map.size() << endl;
   
 
}


void fuse::copyDatasetCtoR( BMAcc_RW<double> &Dnew, const double *Dold ){
   
  uint_fast32_t i,N,K;
  
  double total;
  
  N = Dnew.nrow();
  K = N*N;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,total) shared(Dold,Dnew)
  for(i=0; i<K; i++){
    total   = 0;
    total   = Dold[i];
    Dnew[i] = total;
  }

}



void fuse::transpose( const double *Mat, double *tMat, const vector<tripleInt> &map, const uint_fast32_t n){

  uint_fast32_t i,j,k,kk,N,K,rind;

  N = n;
  K = map.size();
  
#pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(i,j,k,kk,rind) shared(Mat,tMat,map)
  for( kk=0; kk<K; kk++ ){
    
    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
        
    rind = Trag_eq(j,i,N);

    tMat[rind] = 0;
    tMat[rind] = Mat[k];
      
  }   

  
}

void fuse::transpose( const arma::Mat<double> &Mat, arma::Mat<double> &tMat, const uint_fast32_t n, bool upperTri ){
 
  uint_fast32_t i,j,N;
  
  N = n;

  if( upperTri ){

    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j) shared(Mat,tMat)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	tMat.at(j,i) = 0;//Mat.at(i,j);
	tMat.at(i,j) = Mat.at(j,i);
	
      }
      
    }
  }
    
  } else {
    
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j) shared(Mat,tMat)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	tMat.at(j,i) = Mat.at(i,j);
	tMat.at(i,j) = Mat.at(j,i);
	
      }
      
    }
  }
  
  }//ifelse
  
}

void fuse::transpose( const double *Mat, double *tMat, const uint_fast32_t n, bool upperTri ){
 
  uint_fast32_t i,j,N,rind,cind;
  
  N = n;

  if( upperTri ){

    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind) shared(Mat,tMat)
    for(i=0; i<N; i++){
      for(j=0; j<N; j++){    

	if( i <= j ){

	  rind    = (i*N)+j;
	  cind    = (j*N)+i;	
	
	  tMat[cind] = 0;//Mat[rind];
	  tMat[rind] = Mat[cind]; 

	}
      
      }
    }

  } else {
  
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind) shared(Mat,tMat)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	rind    = (i*N)+j;
	cind    = (j*N)+i;	
	
	tMat[cind] = Mat[rind];
	tMat[rind] = Mat[cind]; 

      }
      
    }
  }

  }//ifelse
  
}

void fuse::transpose( BMAcc_RW<double> &Mat, double *tMat, const uint_fast32_t n, bool upperTri ){
 
  uint_fast32_t i,j,N,rind,cind;
  
  N = n;

  if( upperTri ){

    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind) shared(Mat,tMat)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	rind = (i*N)+j;
	cind = (j*N)+i;
	
	tMat[cind] = 0;//Mat[rind];
	tMat[rind] = Mat[cind]; 

      }
      
    }
  }

  } else {
  
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind) shared(Mat,tMat)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	rind = (i*N)+j;
	cind = (j*N)+i;
	
	tMat[cind] = Mat[rind];
	tMat[rind] = Mat[cind]; 

      }
      
    }
  }

  }//ifelse
  
}

void fuse::symmetrise( arma::Mat<double> &Mat, const uint_fast32_t n, double &sum ){
 
  uint_fast32_t i,j,N;
  
  double total;
  
  N = n;

  double *x = (double*)malloc((N*N)*sizeof(double));
  
  //arma::Mat<double> tMat(N,N);
  arma::Mat<double> tMat(x,N,N,false,false);
  
  transpose( Mat, tMat, N, true );

  sum = 0;
  
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,total) shared(Mat,tMat) reduction(+:sum)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){    

      if( i <= j ){

	if( i == j ){
	  total       = 0;
	  total       = Mat.at(i,j)+1;
	  Mat.at(i,j) = total;
	  sum        += total;
	} else {	
	  total       = 0;
	  total       = (tMat.at(i,j) + Mat.at(i,j))/2;
	  Mat.at(i,j) = total;
	  Mat.at(j,i) = total;
	  sum        += 2*total; 
	
	}     
      }

    }
  }
      

  //free space
  if(x){free(x);}
  //reset the number of elements in tMat to zero, space should be freed when leaving function
  //tMat.reset(); 
  
}



void fuse::symmetrise( double *Mat, const uint_fast32_t n, double &sum ){
 
  uint_fast32_t i,j,N,K,rind,cind;
  
  double total;
  
  N = n;
  K = N*N;

  double *tMat = (double*)malloc((N*N)*sizeof(double));
  //double *tMat = (double*)_mm_malloc((N*N)*sizeof(double),CSIZE);  

  transpose( Mat, tMat, N, true );

  sum = 0;

#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind,total) shared(Mat,tMat) reduction(+:sum)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){

      if( i <= j ){

	rind  = (i*N)+j;
	cind  = (j*N)+i;
	total = (Mat[rind]+tMat[rind])/2;
	
	if( i == j ){
	  Mat[rind] = total;
	  sum      += total;
	} else {
	  Mat[rind] = total;
	  Mat[cind] = total;
	  sum      += 2*total;
	}

      }

    }
  }
	

  /*
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,total) shared(Mat,tMat) reduction(+:sum)
  for(i=0; i<K; i++){
    total  = 0;
    total  = (Mat[i]+tMat[i])/2;
    Mat[i] = total;
    sum   += total;
  }
  */
  
  //free space
  if(tMat){free(tMat);}
  //if(tMat){_mm_free(tMat)};
  
  
}

void fuse::symmetrise( BMAcc_RW<double> &Pcnew, double *Mat, const uint_fast32_t n, double &sum ){
 
  uint_fast32_t i,j,N,K,rind,cind;
  
  double total;
  
  N = n;
  K = N*N;
  
  double *tMat = (double*)malloc((N*N)*sizeof(double));
  //double *tMat = (double*)_mm_malloc((N*N)*sizeof(double),CSIZE);  

  transpose( Mat, tMat, N, true);

  sum = 0;

#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind,cind,total) shared(Mat,tMat,Pcnew) reduction(+:sum)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){

      if( i <= j ){

	rind  = (i*N)+j;
	cind  = (j*N)+i;
	total = (Mat[rind]+tMat[rind])/2;
	
	if( i == j ){
	  Pcnew[rind] = total;
	  sum        += total;
	} else {
	  Pcnew[rind] = total;
	  Pcnew[cind] = total;
	  sum        += 2*total;
	}

      }

    }
  }

  /*
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,total) shared(Mat,tMat,Pcnew) reduction(+:sum)
  for(i=0; i<K; i++){
    total    = 0;
    total    = (Mat[i]+tMat[i])/2;
    Pcnew[i] = total;
    sum     += total;
  }
  */
  
  //free space
  if(tMat){free(tMat);}
  //if(tMat){_mm_free(tMat);}
    
}


void fuse::symmetrise( BMAcc_RW<double> &Mnew, const double *Mat, const double *tMat, const uint_fast32_t n ){

  uint_fast32_t i,N,K;
  
  N = n;   
  K = N*N;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i) shared(Mat,tMat,Mnew)
  for(i=0; i<K; i++){
    Mnew[i] = (Mat[i]+tMat[i])/2;
  }

  
}


void fuse::symmetrise( double *Mnew, const double *Mat, const double *tMat, const uint_fast32_t n ){
 
  uint_fast32_t i,N,K;
  
  N = n;   
  K = N*N;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i) shared(Mat,tMat,Mnew)
  for(i=0; i<K; i++){
    Mnew[i] = (Mat[i]+tMat[i])/2;
  }

  
}

void fuse::rowNorm( double *Mnew, const uint_fast32_t n, bool addOne ){
   
  uint_fast32_t i,j,N,rind,cind;
  
  double total;
  
  N = n;

  //row mean
  double *rMn = (double*)malloc(N*sizeof(double));
  //double *rMn = (double*)_mm_malloc(N*sizeof(double),CSIZE);

  //#pragma omp parallel default(none) firstprivate(N) private(i,j,ind) shared(Mnew,rMn)
  //{
    
  //#pragma omp for
    for(i=0; i<N; i++){ rMn[i] = 0; }
  
    #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,rind) shared(Mnew) reduction(+:rMn[:N])
    //#pragma omp for collapse(2) reduction(+:rMn[:N])
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      rind     = (i*N)+j; //row-major
      //cind     = (j*N)+i; //column-major      
      rMn[i] += Mnew[rind];     
    }
  }

 
  
    #pragma omp parallel for schedule(static) default(none) firstprivate(N) private(i,j,rind) shared(Mnew,rMn)
    //#pragma omp for 
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      rind = (i*N)+j;
      //cind = (j*N)+i;
      if( rMn[i] == 0 ){
	Mnew[rind] = 0;
      } else {
	Mnew[rind] /= rMn[i];
      }
    }
  }

  //}//parallel
  
  //free space
  if(rMn){free(rMn);}
  //if(rMn){_mm_free(rMn);}
  
}


void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1, const double *P2, const double *P3, const vector<tupleSr> &map, const vector<tripleInt> &mapP, double &sum ){

  uint_fast32_t i,j,kk,k,KK,KP,K,N,ind,rind;

  double total, SrV, tSrV;

  sum = 0;
  
  N  = n;
  K  = map.size();
  KP = mapP.size();
  KK = N*(N+1)/2;

  //#pragma omp parallel default(none) firstprivate(KK,K,KP) private(i,j,k,kk,SrV,tSrV,total) shared(P1,P2,P3,Pv,map,mapP,sum)
  //{
  
    #pragma omp parallel for schedule(static) default(none) firstprivate(KK) private(k) shared(Pv)
    //#pragma omp for
    for(k=0; k<KK; k++){
    Pv[k] = 0;
  }
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,SrV,tSrV,total) shared(P1,P2,P3,Pv,map)
    //#pragma omp for
  for(kk=0; kk<K; kk++){

    i    = std::get<0>(map[kk]);
    j    = std::get<1>(map[kk]);
    k    = std::get<2>(map[kk]);
    SrV  = std::get<3>(map[kk]);
    tSrV = std::get<4>(map[kk]);   
    
    total = 0;
    total = (P1[k] + P2[k] + P3[k])/3;
    total = tSrV * total * SrV;
    Pv[k] = total; 	

  }
  
  #pragma omp parallel for schedule(static) default(none) firstprivate(KP) private(i,j,k,kk,total) shared(Pv,mapP) reduction(+:sum)
  //#pragma omp for reduction(+:sum)
  for(kk=0; kk<KP; kk++){

    i = std::get<0>(mapP[kk]);
    j = std::get<1>(mapP[kk]);
    k = std::get<2>(mapP[kk]);    

    total   = 0;

    if( i == j ){
      total  = 1;
      Pv[k]  = total;
      sum   += total;
    } else {
      total = Pv[k];
      sum  += 2*total;
    }
    
  }
   
  //}//parallel
  
}

void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1, const vector<tupleSr> &map ){

  uint_fast32_t i,j,kk,k,KK,K,N,ind,rind;

  double total, SrV, tSrV;

  N  = n;
  K  = map.size();
  KK = N*(N+1)/2;


  //#pragma omp parallel default(none) firstprivate(K,KK) private(k,kk,SrV,tSrV,total) shared(P1,P2,P3,Pv,map)
  //{
  
    #pragma omp parallel for schedule(static) default(none) firstprivate(KK) private(k) shared(Pv)
    //#pragma omp for
    for(k=0; k<KK; k++){
    Pv[k] = 0;
  }
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,SrV,tSrV,total) shared(P1,Pv,map)
    //#pragma omp for
  for(kk=0; kk<K; kk++){

    i    = std::get<0>(map[kk]);
    j    = std::get<1>(map[kk]);
    k    = std::get<2>(map[kk]);
    SrV  = std::get<3>(map[kk]);
    tSrV = std::get<4>(map[kk]);   
    
    total = 0;
    total = P1[k];
    total = tSrV * total * SrV;
    Pv[k] = total; 	

  }
  
  //}//parallel
  
}



void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1, const double *P2, const double *P3, const vector<tupleSr> &map ){

  uint_fast32_t i,j,kk,k,KK,K,N,ind,rind;

  double total, SrV, tSrV;

  N  = n;
  K  = map.size();
  KK = N*(N+1)/2;


  //#pragma omp parallel default(none) firstprivate(K,KK) private(k,kk,SrV,tSrV,total) shared(P1,P2,P3,Pv,map)
  //{
  
    #pragma omp parallel for schedule(static) default(none) firstprivate(KK) private(k) shared(Pv)
    //#pragma omp for
    for(k=0; k<KK; k++){
    Pv[k] = 0;
  }
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,SrV,tSrV,total) shared(P1,P2,P3,Pv,map)
    //#pragma omp for
  for(kk=0; kk<K; kk++){

    i    = std::get<0>(map[kk]);
    j    = std::get<1>(map[kk]);
    k    = std::get<2>(map[kk]);
    SrV  = std::get<3>(map[kk]);
    tSrV = std::get<4>(map[kk]);   
    
    total = 0;
    total = (P1[k] + P2[k] + P3[k])/3;
    total = tSrV * total * SrV;
    Pv[k] = total; 	

  }
  
  //}//parallel
  
}

//new code
void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1,
		    CSF *Srcsc, CSF *tSrcsr, const vector<tripleInt> &map, int nTHREADS ){

  
  uint_fast32_t i,j,k,kk,K,N,rind,cind;

  double total, sum;
  
  N = n;
  K = map.size();

  cout << "flg1" << endl;
  
  double *x = (double*)malloc((N*N)*sizeof(double));  

  cout << "flg2" << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(x,P1,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    total  = P1[k];
    
    x[rind] = total;
    x[cind] = total;

  }

  cout << "flg3" << endl;
  
  arma::mat pv(x,N,N,false,false);//no copy of data from x, life span of pv bound to x. 
				      
  cout << "flg4" << endl;

  print_mem_usage();
  
  // res = tSr * pv * Sr
  arma::mat prod, res;
  csr_dense_tcrossprod( tSrcsr, pv, prod, nTHREADS );

  cout << "flg5" << endl;

  print_mem_usage();    
  
  //free space, this will free pv as well.
  if(x){free(x); }

  cout << "flg6" << endl; 
  
  print_mem_usage();    
  
  dense_csc_prod      ( prod, Srcsc, res, nTHREADS );

  cout << "flg7" << endl; 

  print_mem_usage();
    
   //arma::mat objects should be deleted once out of scope of function, but use reset function
  //to set elements in these matrices to zero.
  prod.reset();  

  cout << "flg8" << endl;

  cout << "is res symmetric: " << res.is_symmetric() << endl;

  if( !res.is_symmetric() ){
    cout << " (Normalising... ";      
    symmetrise( res, N, sum );
    cout << " (" << sum << ") done.)" << endl;
  }

  cout << "is res symmetric: " << res.is_symmetric() << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(Pv,res,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    Pv[k] = 0;
    Pv[k] = res.at(i,j);
    
  }

  //reset res
  res.reset();
  
  cout << "flg9" << endl;
   
  
}


//new code
void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1, const double *P2, const double *P3, CSF *Srcsc, CSF *tSrcsr, const vector<tripleInt> &map, int nTHREADS ){

  
  uint_fast32_t i,j,k,kk,K,N,rind,cind;

  double total, sum;
  
  N = n;
  K = map.size();

  cout << "flg1" << endl;
  
  double *x = (double*)malloc((N*N)*sizeof(double));  

  cout << "flg2" << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(x,P1,P2,P3,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    total  = (P1[k]+P2[k]+P3[k])/3;
    
    x[rind] = total;
    x[cind] = total;

  }

  cout << "flg3" << endl;
  
  arma::mat pv(x,N,N,false,false);//no copy of data from x, life span of pv bound to x. 
				      
  cout << "flg4" << endl;

  print_mem_usage();
  
  // res = tSr * pv * Sr
  arma::mat prod, res;
  csr_dense_tcrossprod( tSrcsr, pv, prod, nTHREADS );

  cout << "flg5" << endl;

  print_mem_usage();    
  
  //free space, this will free pv as well.
  if(x){free(x); }

  cout << "flg6" << endl; 
  
  print_mem_usage();    
  
  dense_csc_prod      ( prod, Srcsc, res, nTHREADS );

  cout << "flg7" << endl; 

  print_mem_usage();
    
   //arma::mat objects should be deleted once out of scope of function, but use reset function
  //to set elements in these matrices to zero.
  prod.reset();  

  cout << "flg8" << endl;

  cout << "is res symmetric: " << res.is_symmetric() << endl;

  if( !res.is_symmetric() ){
    cout << " (Normalising... ";      
    symmetrise( res, N, sum );
    cout << " (" << sum << ") done.)" << endl;
  }

  cout << "is res symmetric: " << res.is_symmetric() << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(Pv,res,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    Pv[k] = 0;
    Pv[k] = res.at(i,j);
    
  }

  //reset res
  res.reset();
  
  cout << "flg9" << endl;
   
  
}


//new code
void fuse::updateP( double *Pv, const uint_fast32_t n, const double *P1, const double *P2, const double *P3, const vector<tupleCOO> &Srmap, const vector<tripleInt> &map ){

  
  uint_fast32_t i,j,k,kk,K,N,rind,cind;

  double total;
  
  N = n;
  K = map.size();

  cout << "flg1" << endl;
  
  double *x = (double*)malloc((N*N)*sizeof(double));  

  cout << "flg2" << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(x,P1,P2,P3,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    total  = (P1[k]+P2[k]+P3[k])/3;
    
    x[rind] = total;
    x[cind] = total;

  }

  cout << "flg3" << endl;
  
  arma::mat pv(x,N,N,false,false);//no copy of data from x, life span of pv bound to x. 
				      
  cout << "flg4" << endl;
				      
  // res = tSr * pv * Sr
  arma::mat res;

  arma::SpMat<double> Sr; Sr.zeros(N,N);
  SrVec2ArmaD( Sr, Srmap, N );

  cout << "flg5" << endl;
  
  res = Sr.t() * pv * Sr;

  cout << "flg6" << endl;
  
#pragma omp parallel for schedule(static) firstprivate(K) private(i,j,k,kk,rind,cind,total) shared(Pv,res,map)
  for(kk=0; kk<K;kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);
    
    rind = (i*N)+j;
    cind = (j*N)+i;     
	
    Pv[k] = 0;
    Pv[k] = res(i,j);
    
  }

  cout << "flg7" << endl;
  
  //free space
  if(x){free(x); }

  cout << "flg8" << endl;
  
  //arma::mat objects should be deleted once out of scope of function, but use reset function
  //to set elements in these matrices to zero.
  Sr.reset();
  res.reset();

  cout << "flg9" << endl;
  
}

void fuse::convergePc( const double *P1, const double *P2,
		       const vector<tripleInt> &map, double &sum, const uint_fast32_t n ){

  uint_fast32_t i,j,k,kk,N,K;

  double total;
  
  sum = 0;

   N = n;
   K = map.size();
   
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,total) shared(P1,P2,map) reduction(+:sum)
   for(kk=0; kk<K; kk++){

     i = std::get<0>(map[kk]);
     j = std::get<1>(map[kk]);
     k = std::get<2>(map[kk]);
     
     total = 0;
     total = (P1[k] + P2[k])/2;
     
     if( i == j ){
       sum += total;
     } else {
       sum += 2*total;
     }
     
   }
  
}


void fuse::convergePc( const double *P1, const double *P2, const double *P3, const double *P4,
		       const vector<tripleInt> &map, double &sum, const uint_fast32_t n ){

  uint_fast32_t i,j,k,kk,N,K;

  double total;
  
  sum = 0;

   N = n;
   K = map.size();
   
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,total) shared(P1,P2,P3,P4,map) reduction(+:sum)
   for(kk=0; kk<K; kk++){

     i = std::get<0>(map[kk]);
     j = std::get<1>(map[kk]);
     k = std::get<2>(map[kk]);
     
     total = 0;
     total = (P1[k] + P2[k] + P3[k] + P4[k])/4;
     
     if( i == j ){
       sum += total;
     } else {
       sum += 2*total;
     }
     
   }
  
}

 
//this code should be fine, as P matrices are symmetric
void fuse::normaliseP( double *Pnew, const vector<tripleInt> &map, const uint_fast32_t n, double &sum ){

  
  uint_fast32_t i,j,k,kk,N,K;  

  double total;
  
  N = n;
  K = map.size();

  sum = 0;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,k,kk,total) shared(Pnew,map) reduction(+:sum)
  for( kk=0; kk<K; kk++ ){
    
    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);

    total = 0;
    
    if( i == j ){
      total    = 1;
      Pnew[k] += total;
      sum     += total;
    } else {
      total    = Pnew[k];
      sum     += 2*total;
    }
       
      
  }

 
  //symmetrise(Pnew,map,N,sum);

 }


void fuse::normaliseP( double *Pnew, const uint_fast32_t n, double &sum ){
 
  uint_fast32_t i,j,N,K,ind,rind;
  
  double total; 
  
  N = n;
  K = N*N;
  
  #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,ind) shared(Pnew)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){

      if( i <= j ){

	ind = (i*N)+j;

	if( i == j ){
	  Pnew[ind] += 1;	 	  
	}
	 
	
      }
    }
  }  


  //symmetrise
  symmetrise( Pnew, N, sum );  

}

void fuse::copyP( const double *P1new, const double *P2new,
		  double *P1old, double *P2old,
		  const vector<tripleInt> &map, const uint_fast32_t n){

  uint_fast32_t k,K;
  
  K = n;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(k) shared(P1old,P2old,P1new,P2new)
  for(k=0; k<K; k++){
	P1old[k] = P1new[k];
	P2old[k] = P2new[k];
  }

}


void fuse::copyP( const double *P1new, const double *P2new, const double *P3new, const double *P4new,
		  double *P1old, double *P2old, double *P3old, double *P4old,
		  const vector<tripleInt> &map, const uint_fast32_t n){

  uint_fast32_t k,K;
  
  K = n;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(k) shared(P1old,P2old,P3old,P4old,P1new,P2new,P3new,P4new)
  for(k=0; k<K; k++){
	P1old[k] = P1new[k];
	P2old[k] = P2new[k];
	P3old[k] = P3new[k];
	P4old[k] = P4new[k];
  }

}


void fuse::copyP( const double *P1new, const double *P2new, const double *P3new, const double *P4new,
	    double *P1old, double *P2old, double *P3old, double *P4old, const uint_fast32_t n){

 
  uint_fast32_t i,N,K;
  
  N = n;   
  K = N*N;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i) shared(P1old,P2old,P3old,P4old,P1new,P2new,P3new,P4new)
  for(i=0; i<K; i++){
	P1old[i] = P1new[i];
	P2old[i] = P2new[i];
	P3old[i] = P3new[i];
	P4old[i] = P4new[i];
  }

}

void fuse::updatePc( BMAcc_RW<double> &Pc, const vector<tripleInt> &map,
		     const double *P1, const double *P2, double &sum ){

  uint_fast32_t i,j,k,kk,N,K,rind,cind;
  
  double total;
  
  N = Pc.nrow();
  K = map.size();
  
  sum = 0;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K,N) private(i,j,k,kk,rind,cind,total) shared(P1,P2,Pc,map) reduction(+:sum)
  for(kk=0; kk<K; kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);

    rind = (i*N)+j;
    cind = (j*N)+i;

    total    = (P1[k] + P2[k])/2;    
    Pc[rind] = total;
    Pc[cind] = total;

    if( i == j ){    
      sum   += total;
    } else {
      sum   += 2*total;
    }
    
  }
   
  
}


void fuse::updatePc( BMAcc_RW<double> &Pc, const vector<tripleInt> &map, const double *P1, const double *P2,
	       const double *P3, const double *P4, double &sum ){

  uint_fast32_t i,j,k,kk,N,K,rind,cind;
  
  double total;
  
  N = Pc.nrow();
  K = map.size();
  
  sum = 0;
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K,N) private(i,j,k,kk,rind,cind,total) shared(P1,P2,P3,P4,Pc,map) reduction(+:sum)
  for(kk=0; kk<K; kk++){

    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);

    rind = (i*N)+j;
    cind = (j*N)+i;

    total    = (P1[k] + P2[k] + P3[k] + P4[k])/4;    
    Pc[rind] = total;
    Pc[cind] = total;

    if( i == j ){    
      sum   += total;
    } else {
      sum   += 2*total;
    }
    
  }
   
  
}



void fuse::normalisePc( BMAcc_RW<double> &Pnew, const uint_fast32_t n, double &sum ){

 
  uint_fast32_t i,j,N,K,ind,rind;
  
  double total;

  //sum = 0;
  
  N = n;
  K = N*N;

  double *Ptmp = (double*)malloc((N*N)*sizeof(double));
  //double *Ptmp = (double*)_mm_malloc((N*N)*sizeof(double),CSIZE);  
  
  //copy R matrix to C array
  copyDatasetRtoC( Ptmp, Pnew );    
  
  //old sequence (start)
  // row normalisation
  rowNorm( Ptmp, N ); 

  //symmetrise
  symmetrise( Pnew, Ptmp, N, sum );
  //old sequence (end)

  //free space
  if(Ptmp){free(Ptmp);}
  //if(Ptmp){_mm_free(Ptmp);}
  
}

void fuse::map2Array( const double *Mat, double *Arr, const uint_fast32_t n, const vector<tripleInt> &map ){

  uint_fast32_t i,j,k,kk,ind,N,K;

  N = n;
  K = map.size();
  
#pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(i,j,k,kk,ind) shared(Mat,Arr,map)
  for( kk=0; kk<K; kk++ ){
    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);

    ind = (i*N)+j;

    Arr[k] = Mat[ind];
    
  }   
  
}


//indcies which map from 2d matrix (with Diagonal) to 1d array
void fuse::mapMatrixDiag2Array( const uint_fast32_t N, vector<tripleInt> &map ){

  uint_fast32_t i,j,k,K;

  if( N > 0 ){

    K = map.size();   


    // user-defined reduction for pushing back in vector map
    #pragma omp declare reduction (merge : std::vector<tripleInt> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N) private(i,j,k) reduction(merge:map)
  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){

      if( i <= j ){
      
	k = Trag_eq(i,j,N);

	if( k != -1 ){

	  map.push_back( tripleInt(i,j,k) );
	  
	  //std::get<0>(map[k]) = i;
	  //std::get<1>(map[k]) = j;
	  //std::get<2>(map[k]) = k;

	  //} else {

	  //std::get<0>(map[k]) = -1;
	  //std::get<1>(map[k]) = -1;
	  //std::get<2>(map[k]) = -1;

	  //}

	}
      }
    }
  }
 

  }
  
}

//indcies which map from 1d array to 2d matrix (with Diagonal)
void fuse::mapArray2MatrixDiag( const uint_fast32_t N, vector<tripleInt> &map ){

  uint_fast32_t i,j,k,K;

  if( N > 0 ){

    K = map.size();   

    //init the map
    for( k=0; k<K; k++ ){
      map.push_back( tripleInt(0,0,k) );
    }
  
#pragma omp parallel for schedule(static) default(none) firstprivate(N,K) private(i,j,k) shared(map)
  for( k=0; k<K; k++ ){    

    i = 0;
    j = 0;
    
    Trag_reverse_eq(k, N, i, j);
    
    std::get<0>(map[k]) = i;
    std::get<1>(map[k]) = j;   
    
  }
     

  }
  
  
}

