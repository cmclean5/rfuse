#include "Headers.h"
#include "fuse.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd rblock/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
// If need to run (to build Makevars file in src):
// 1) $> cd rblock/
// 2) $> autoconfig
// 3) $> ./configure
// BUILD PACKAGE:
// 1) $> R CMD build rblock/
// INSTALL PACKAGE (note 'INSTALL' should be in capital letters):
// 1) $> R CMD INSTALL rblock_1.0.tar.gz
//----------------

//Global;

fuse *model=nullptr;
double OMP_start;
double OMP_end;

//------ declare R functions ------

// [[Rcpp::export]]
void createModel(){ model = new fuse(); }

// [[Rcpp::export]]
void deleteModel(){ if(model!=nullptr){ delete model;} }

// [[Rcpp::export]]
void freeCspace(){if(model!=nullptr){ model->freeSpace(); } }

  
// [[Rcpp::export]]
void startTiming(){

  if(model){ cout << "create model first." << endl; }
  else {  
  model->getOpenMPStart(OMP_start);
  }
}

// [[Rcpp::export]]
void endTiming( IntegerVector print=1 ){  

  if(model==nullptr){ cout << "create model first." << endl; return; }

  model->getOpenMPStart(OMP_end);

  if( print[0] == 1 ){
    cout << "Work took " << (OMP_end-OMP_start) << " seconds" << endl;
  }
  
}


// [[Rcpp::export]]
Rcpp::List getEdgelist2( Environment Sfuse, Environment Mat, IntegerVector nCORES, IntegerVector Norm2=0 ){
 
  // read only the matrix 
  XPtr<FBM> xSc = Sfuse["address"];
  BMAcc<double> Sc(xSc);

  // read only the matrix 
  XPtr<FBM> xMat = Mat["address"];
  BMAcc<double> M(xMat);  

  Rcpp::NumericMatrix ed;
  
  if(model==nullptr){
    cout << "create model first." << endl;
    ed (0,3);
    return Rcpp::List::create(Rcpp::Named("RESULT")=ed);
  }
    
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  uint_fast32_t i,j,k,kk,ind,N,K,Ks,norm2;
  double cand_edge, weight;
  bool useLog=false;
  
  N     = M.nrow();  
  //K     = N*(N+1)/2; //with    diagonal
  K     = N*(N-1)/2; //without diagonal
  norm2 = Norm2[0];  
  Ks    = 0;
  
  vector<EDGE> edges;
  edges.reserve(K);  

  cout << "edges capacity: " << edges.capacity() << " edges.size() " << edges.size() << endl;
  

  // user-defined reduction for pushing back in vector map
  #pragma omp declare reduction (merge : std::vector<EDGE> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

   
  #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N,norm2) private(i,j,ind,cand_edge, weight) shared(M,Sc) reduction(merge: edges)
  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){

      if( i < j ){

	  ind       = (i*N)+j;	
	  cand_edge = M[ind];
	  weight    = Sc[ind];
	
	  if( norm2 == 1 ){ weight /= 0.5; }
	
	  if( cand_edge != 0 ){
	    edges.push_back( EDGE((i+1),(j+1), weight) );	    
	  } 
      }
    
    }
  }     

  cout << "edges capacity: " << edges.capacity() << " edges.size() " << edges.size() << endl;
  
  ed = Rcpp::NumericMatrix(edges.size(),3);
  for( k=0; k<edges.size(); k++ ){
      ed(k,0) = std::get<0>(edges[k]);
      ed(k,1) = std::get<1>(edges[k]);
      ed(k,2) = std::get<2>(edges[k]);
  }

  //free edges memory vector
  std::vector<EDGE>().swap(edges);
  
  return Rcpp::List::create(Rcpp::Named("RESULT")=ed);


}


// [[Rcpp::export]]
Rcpp::List getEdgelist( Environment Min, IntegerVector nCORES, IntegerVector Norm2=0 ){
 
  // read only the matrix 
  XPtr<FBM> xMin = Min["address"];
  BMAcc<double> M(xMin);  

  Rcpp::NumericMatrix ed;
  
  if(model==nullptr){
    cout << "create model first." << endl;
    ed (0,3);
    return Rcpp::List::create(Rcpp::Named("RESULT")=ed);
  }
    
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  uint_fast32_t i,j,k,kk,ind,N,K,Ks,norm2;
  double weight;
  bool useLog=false;
  
  N     = M.nrow();  
  //K     = N*(N+1)/2; //with    diagonal
  K     = N*(N-1)/2; //without diagonal
  norm2 = Norm2[0];  
  Ks    = 0;
  
  vector<EDGE> edges;
  edges.reserve(K);  

  cout << "edges capacity: " << edges.capacity() << " edges.size() " << edges.size() << endl;
  

  // user-defined reduction for pushing back in vector map
  #pragma omp declare reduction (merge : std::vector<EDGE> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

   
#pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(N,norm2) private(i,j,ind,weight) shared(M) reduction(merge: edges)
  for( i=0; i<N; i++ ){
    for( j=0; j<N; j++ ){

      if( i < j ){

	  ind     = (i*N)+j;	
	  weight  = M[ind];
	
	  if( norm2 == 1 ){ weight /= 0.5; }
	
	  if( weight != 0 ){
	    edges.push_back( EDGE((i+1),(j+1), weight) );	    
	  } 
      }
    
    }
  }     

  cout << "edges capacity: " << edges.capacity() << " edges.size() " << edges.size() << endl;

  ed = Rcpp::NumericMatrix(edges.size(),3);
  for( k=0; k<edges.size(); k++ ){
      ed(k,0) = std::get<0>(edges[k]);
      ed(k,1) = std::get<1>(edges[k]);
      ed(k,2) = std::get<2>(edges[k]);
  }

  //free edges memory vector
  std::vector<EDGE>().swap(edges);
  
  return Rcpp::List::create(Rcpp::Named("RESULT")=ed);


}

// [[Rcpp::export]]
void tcpp(Environment Min, NumericVector nCORES){

  // pointer to matrix object on disk
  XPtr<FBM_RW> xMin = Min["address_rw"];
  BMAcc_RW<double> M(xMin);

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );
  
  uint_fast32_t N = M.nrow();  
  
  double *Mtmp  = (double*)malloc((N*N)*sizeof(double));
  double *tMtmp = (double*)malloc((N*N)*sizeof(double)); 
  

  cout << "M ";  model->copyDatasetRtoC(Mtmp, M);          cout << "...done." << endl;
  cout << "tM "; model->transpose(Mtmp, tMtmp, N, false);  cout << "...done." << endl;

  model->symmetrise( M, Mtmp, tMtmp, N );

  //free space
  if(Mtmp){free(Mtmp);}
  if(tMtmp){free(tMtmp);}

 
}



// [[Rcpp::export]]
void readPmatrix(Environment Pmat, NumericVector Dataset, NumericVector nCORES){
  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xPmat = Pmat["address_rw"];
  BMAcc_RW<double> Pr(xPmat);

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  //loop index
  uint_fast32_t k;
  
  //size of Smat
  uint_fast32_t N = Pr.nrow();  

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K);  
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }
  
  double sum = 0.0;  
  
  
  try{
    
    switch (dataset) {
    case 1:
      cout << "P_bms ";    
      model->PrDiag1 = (double*)malloc(K*sizeof(double));
      model->copyDatasetRtoC( Pr, model->PrDiag1, N, model->map2Arr, sum ); cout << sum << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 2:
      cout << "P_scq ";
      model->PrDiag2 = (double*)malloc(K*sizeof(double));
      model->copyDatasetRtoC( Pr, model->PrDiag2, N, model->map2Arr, sum ); cout << sum << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 3:
      cout << "P_rbsr ";
      model->PrDiag3 = (double*)malloc(K*sizeof(double));
      model->copyDatasetRtoC( Pr, model->PrDiag3, N, model->map2Arr, sum ); cout << sum << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 4:
      cout << "P_dcdq ";
      model->PrDiag4 = (double*)malloc(K*sizeof(double));
      model->copyDatasetRtoC( Pr, model->PrDiag4, N, model->map2Arr, sum ); cout << sum << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;

    default:
      cout << "parameter Dataset needs to be in range [1,4]." << endl;

    }

       
      
  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;    
  }  

    
  //free mapping vector
  //std::vector<tripleInt>().swap(map2Arr);
  //freeMapSpace();
  
  
}


// [[Rcpp::export]]
void readSmatrix2(Environment Smat, NumericVector Dataset, NumericVector nCORES){
  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xSmat = Smat["address_rw"];
  BMAcc_RW<double> Sr(xSmat);

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  uint_fast32_t k;
  
  //size of Smat
  uint_fast32_t N = Sr.nrow();  

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K);  
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }

  
  try{
    
    switch (dataset) {
    case 1:
      cout << "Sr_bms ";
      model->tSr1csr = new CSF();
      model->Sr1csc  = new CSF();
      model->mapSrmatrix( Sr, N, model->tSr1csr, model->Sr1csc, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 2:
      cout << "Sr_scq ";    
      model->tSr2csr = new CSF();
      model->Sr2csc  = new CSF();
       model->mapSrmatrix( Sr, N, model->tSr2csr, model->Sr2csc, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 3:
      cout << "Sr_rbsr ";
      model->tSr3csr = new CSF();
      model->Sr3csc  = new CSF();
      model->mapSrmatrix( Sr, N, model->tSr3csr, model->Sr3csc, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 4:
      cout << "Sr_dcdq ";    
      model->tSr4csr = new CSF();
      model->Sr4csc  = new CSF();
      model->mapSrmatrix( Sr, N, model->tSr4csr, model->Sr4csc, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;

    default:
      cout << "Dataset needs to be in range [1,4]." << endl;
    }
    
    
  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;    
  }

  
  //free mapping vector
  //std::vector<tripleInt>().swap(map2Arr);
  //freeMapSpace();
  
}




// [[Rcpp::export]]
void readSmatrix(Environment Smat, NumericVector Dataset, NumericVector nCORES){

  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xSmat = Smat["address_rw"];
  BMAcc_RW<double> Sr(xSmat);

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  uint_fast32_t k;
  
  //size of Smat
  uint_fast32_t N = Sr.nrow();  

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }

  
  try{
    
    switch (dataset) {
    case 1:
      cout << "Sr_bms ";     
      model->mapSrmatrix( Sr, N, model->mapSr1, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 2:
      cout << "Sr_scq ";    
      model->mapSrmatrix( Sr, N, model->mapSr2, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 3:
      cout << "Sr_rbsr ";    
      model->mapSrmatrix( Sr, N, model->mapSr3, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 4:
      cout << "Sr_dcdq ";    
      model->mapSrmatrix( Sr, N, model->mapSr4, model->map2Arr );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;

    default:
      cout << "Dataset needs to be in range [1,4]." << endl;
    }
    
    
  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;    
  }

  
  //free mapping vector
  //std::vector<tripleInt>().swap(map2Arr);
  //freeMapSpace();
  
}

/*
//// [[Rcpp::export]]
void readSmatrix3(Environment Smat, NumericVector Dataset, NumericVector nCORES){

  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xSmat = Smat["address_rw"];
  BMAcc_RW<double> Sr(xSmat);

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  uint_fast32_t k;
  
  //size of Smat
  uint_fast32_t N = Sr.nrow();  

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }

  
  try{
    
    switch (dataset) {
    case 1:
      cout << "Sr_bms ";     
      model->mapSrmatrix( Sr, N, model->COOSr1, model->map2Arr, true );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 2:
      cout << "Sr_scq ";    
      model->mapSrmatrix( Sr, N, model->COOSr2, model->map2Arr, true );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 3:
      cout << "Sr_rbsr ";    
      model->mapSrmatrix( Sr, N, model->COOSr3, model->map2Arr, true );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;
    case 4:
      cout << "Sr_dcdq ";    
      model->mapSrmatrix( Sr, N, model->COOSr4, model->map2Arr, true );
      cout << " ...done." << endl;
      Rcpp::checkUserInterrupt();
      break;

    default:
      cout << "Dataset needs to be in range [1,4]." << endl;
    }
    
    
  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;    
  }

  
  //free mapping vector
  //std::vector<tripleInt>().swap(map2Arr);
  //freeMapSpace();
  
}
*/

// [[Rcpp::export]]
void normalisePc(Environment Pcmat, NumericVector nCORES){
  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xPcmat = Pcmat["address_rw"];
  BMAcc_RW<double> Pc(xPcmat);
  
  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );  
  
  //size of Smat
  uint_fast32_t N = Pc.nrow();  

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  double sum = 0.0;
  
  try{

    // ensure affinity matrix is normalised and symmetrical
    cout << "run normalisePc: " << endl;
    cout << "Pc "; model->normalisePc( Pc, N, sum ); cout << sum << " ...done." << endl;      
    
 } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;    
  }
  
  
  
}


// [[Rcpp::export]]
void test( Environment Pmat, Environment Smat, NumericVector Dataset, NumericVector nCORES){
  
  
  // pointer to matrix object on disk
  XPtr<FBM_RW> xSmat = Smat["address_rw"];
  BMAcc_RW<double> Sr(xSmat); 

  if(model==nullptr){ cout << "create model first." << endl; return; }
  
  arma::mat Y  = FBM2arma( Pmat ); 

  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  uint_fast32_t i,j,k;

  double val;
  
  //size of Smat
  uint_fast32_t N = Sr.nrow();
  
  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }

  //calculate t(Sr) * Y * Sr => t(CSFcsr) * t(Y) = res, res2 = res * CSFcsc  
  vector<tupleCOO> Srmap, tSrmap;
 
  CSF* Srcsr  = new CSF(); //store Sr1 transpose in csr format
  CSF* Srcsc  = new CSF(); //store Sr1           in csc format
  
  model->mapSrmatrix( Sr, N, Srmap, model->map2Arr, true );//read Sr in column-major format

  tSrVec( Srmap, tSrmap, N );//transpose

  cout << "Y" << endl;
  cout << Y << endl;

  cout << "tSrmap" << endl;
  for(k=0; k<tSrmap.size(); k++){
    i = std::get<0>(tSrmap[k]);
    j = std::get<1>(tSrmap[k]);
    val = std::get<2>(tSrmap[k]);
    cout << "(" << i << "," << j << ") = " << val << endl;
  }

    
  model->mappingCOO2CSF( tSrmap, Srcsr, N, false );//store Sr in CSR format
  model->mappingCOO2CSF(  Srmap, Srcsc, N, true  );//store Sr in CSC format
 
  cout << "Srcsr" << endl;
  for(i=0; i<Srcsr->nz;i++){
    cout << "[" << i << "] val: " << Srcsr->val[i] << ", col: " << Srcsr->col[i] << endl; 
  }

  for(i=0; i<(Srcsr->nrows+1);i++){
    cout << "[" << i << "] row: " << Srcsr->row[i] << endl; 
  }
  cout << "---" << endl;


  arma::mat res,res2;
  model->csr_dense_tcrossprod( Srcsr,   Y, res   );
  model->dense_csc_prod      ( res, Srcsc, res2  );

  cout << "tSr * Y * Sr" << endl;
  cout << res2 << endl;


  arma::SpMat<double> X; X.zeros(N,N);
  SrVec2ArmaD( X, Srmap, N );

  arma::mat res3;
  cout << "tSr * Y * Sr arma " << endl;
  res3 = X.t() * Y * X;
  cout << res3 << endl;

  bool same = approx_equal(res2, res3, "reldiff", 0.1);

  if( same ){
    cout << "matices are the same" << endl;
  } else {
    cout << "matices are different" << endl;
  }
   

  //free space
  model->freeMapSpace();
  std::vector<tupleCOO>().swap(Srmap);
  std::vector<tupleCOO>().swap(tSrmap); 
  if(Srcsr){ delete Srcsr; }
  if(Srcsc){ delete Srcsc; }
 
  
}

// [[Rcpp::export]]  
Rcpp::List findASDedges(Environment Pfuse,IntegerVector nCORES, IntegerVector nSTEPS=0 ){
  //code used to determine number of edges  
    
  // pointer to Pc matrix object on disk
  XPtr<FBM_RW> pc = Pfuse["address_rw"];
  BMAcc_RW<double> Pc(pc);

  Rcpp::NumericMatrix OO;
   
  if(model==nullptr){
    cout << "create model first." << endl;
    OO =  Rcpp::NumericMatrix(0,4);
    return Rcpp::List::create(Rcpp::Named("RESULT")=OO);    
  }

  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] ); 
  
  //---Setup fusion convergous 
  double FUSE_ACC     = 1.0e-6;
  int    FUSE_MAXSTEP = nSTEPS[0]+1;
  int    steps        = 0; 
  double delta        = 0;
  double delta_old    = 0;
  double diff         = 1;
  double sum          = 0;
  bool runPcNorm      = true;  

  std::cout << "Max number of steps: " << FUSE_MAXSTEP << endl;
  
  //set output
  OO = Rcpp::NumericMatrix((FUSE_MAXSTEP+1),4);
  for(int i=0; i<(FUSE_MAXSTEP+1); i++){
    OO(i,0) = 0;
    OO(i,1) = 0;
    OO(i,2) = 0;
    OO(i,3) = 0;
  } 
  
  //set matrix size
  uint_fast32_t N = Pc.nrow(); 

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  if( model->map2Arr.empty() ){
  
    //mapping vector, make sure we reserve correct amount of space.
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 

    //get a mapping from 2d matrix (with diagonal) to 1d array
    model->mapMatrixDiag2Array( N, model->map2Arr ); 

  }
  
  double *P1new = (double*)malloc(K*sizeof(double));
  double *P2new = (double*)malloc(K*sizeof(double));
  double *P3new = (double*)malloc(K*sizeof(double));
  double *P4new = (double*)malloc(K*sizeof(double));
  

  try{
     
      // check convergence 
      cout << "run updatePc... ";
      model->convergePc( model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, sum, N );
      Rcpp::checkUserInterrupt();
      cout << "done." << endl;     
      
      delta = sum;
      diff  = fabs(delta-delta_old);
  
      cout << "> diffusion step: " << 0 << " sum: " << sum << " sum norm: " << sum/(N*N)
	   << endl;

      OO(0,0) = 0;
      OO(0,1) = sum;
      OO(0,2) = sum/(N*N);
      OO(0,3) = diff;

    
      // perform the diffusion for t iterations
      cout << "perform diffusion:" << endl;
      for( steps=1; steps <= FUSE_MAXSTEP; steps++ ){
  
	// update each dataset and normalize each new obtained dataset.
	cout << "run updateP:" << endl;	
	cout << "P1new "; model->updateP( P1new, N, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->mapSr1 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P2new "; model->updateP( P2new, N, model->PrDiag1, model->PrDiag3, model->PrDiag4, model->mapSr2 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P3new "; model->updateP( P3new, N, model->PrDiag1, model->PrDiag2, model->PrDiag4, model->mapSr3 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "P4new "; model->updateP( P4new, N, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->mapSr4 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "done." << endl;

	// check convergence 
	cout << "run convergePc... ";
	model->convergePc( P1new, P2new, P3new, P4new, model->map2Arr, sum, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;     

	delta_old = delta;
	delta     = sum;
	diff      = fabs(delta-delta_old);
	
	cout << "> diffusion step: " << steps << " sum: " << sum << " sum norm: "
	     << sum/(N*N) << " diff: " << diff << endl;
	

	OO(steps,0) = steps;
	OO(steps,1) = sum;
	OO(steps,2) = sum/(N*N);
	OO(steps,3) = diff;


	// Normalize each new obtained networks.
	cout << "run normaliseP: " << endl;
	cout << "P1new "; model->normaliseP( P1new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P2new "; model->normaliseP( P2new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P3new "; model->normaliseP( P3new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P4new "; model->normaliseP( P4new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "done." << endl; 
	
	cout << "run copyP... ";    
	model->copyP( P1new, P2new, P3new, P4new, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;    
	
    
	}//for
  
    

  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;
    runPcNorm = false;
  }
 

  //free memory
  model->freePrDiagSpace();
  model->freeSrSpace();
  
  if( runPcNorm ){
  
    try{
    
      // construct the combined affinity matrix by summing diffused matrices
      cout << "run updatePc... "; 
      model->updatePc( Pc, model->map2Arr, P1new, P2new, P3new, P4new, sum );
      Rcpp::checkUserInterrupt();
      cout << " done." << endl;
      
      // ensure affinity matrix is normalised and symmetrical
      ////cout << "run normalisePc: " << endl;
      ////cout << "Pc "; normalisePc( Pc, N, sum ); cout << sum << " ...done." << endl;
      
    } catch (Rcpp::internal::InterruptedException& e)
      {
	Rcout << "Caught an interrupt!" << std::endl;
      }
    
  }//if

  //free mapping vector
  model->freeMapSpace();
  
  
  //free Pnew space 
  if(P1new){free(P1new);}
  if(P2new){free(P2new);}
  if(P3new){free(P3new);}
  if(P4new){free(P4new);}

  
  //check we've freed space
  model->freeSpace();
   
  
  return Rcpp::List::create(Rcpp::Named("RESULT")=OO);

  
}


// [[Rcpp::export]]  
Rcpp::List findCONTROLedges(Environment Pfuse,IntegerVector nCORES, IntegerVector nSTEPS=0 ){
  //code used to determine number of edges  
    
  // pointer to Pc matrix object on disk
  XPtr<FBM_RW> pc = Pfuse["address_rw"];
  BMAcc_RW<double> Pc(pc);

  Rcpp::NumericMatrix OO;
   
  if(model==nullptr){
    cout << "create model first." << endl;
    OO =  Rcpp::NumericMatrix(0,4);
    return Rcpp::List::create(Rcpp::Named("RESULT")=OO);    
  }

  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] ); 
  
  //---Setup fusion convergous 
  double FUSE_ACC     = 1.0e-6;
  int    FUSE_MAXSTEP = nSTEPS[0]+1;
  int    steps        = 0; 
  double delta        = 0;
  double delta_old    = 0;
  double diff         = 1;
  double sum          = 0;
  bool runPcNorm      = true;  

  std::cout << "Max number of steps: " << FUSE_MAXSTEP << endl;
  
  //set output
  OO = Rcpp::NumericMatrix((FUSE_MAXSTEP+1),4);
  for(int i=0; i<(FUSE_MAXSTEP+1); i++){
    OO(i,0) = 0;
    OO(i,1) = 0;
    OO(i,2) = 0;
    OO(i,3) = 0;
  } 
  
  //set matrix size
  uint_fast32_t N = Pc.nrow(); 

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  if( model->map2Arr.empty() ){
  
    //mapping vector, make sure we reserve correct amount of space.
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 

    //get a mapping from 2d matrix (with diagonal) to 1d array
    model->mapMatrixDiag2Array( N, model->map2Arr ); 

  }
  
  double *P1new = (double*)malloc(K*sizeof(double));
  double *P2new = (double*)malloc(K*sizeof(double));

  try{
     
      // check convergence 
      cout << "run updatePc... ";
      model->convergePc( model->PrDiag1, model->PrDiag2, model->map2Arr, sum, N );
      Rcpp::checkUserInterrupt();
      cout << "done." << endl;     
      
      delta = sum;
      diff  = fabs(delta-delta_old);
  
      cout << "> diffusion step: " << 0 << " sum: " << sum << " sum norm: " << sum/(N*N)
	   << endl;

      OO(0,0) = 0;
      OO(0,1) = sum;
      OO(0,2) = sum/(N*N);
      OO(0,3) = diff;

    
      // perform the diffusion for t iterations
      cout << "perform diffusion:" << endl;
      for( steps=1; steps <= FUSE_MAXSTEP; steps++ ){
  
	// update each dataset and normalize each new obtained dataset.
	cout << "run updateP:" << endl;	
	cout << "P1new "; model->updateP( P1new, N, model->PrDiag2, model->mapSr1 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P2new "; model->updateP( P2new, N, model->PrDiag1, model->mapSr2 ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();      
	cout << "done." << endl;

	// check convergence 
	cout << "run convergePc... ";
	model->convergePc( P1new, P2new, model->map2Arr, sum, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;     

	delta_old = delta;
	delta     = sum;
	diff      = fabs(delta-delta_old);
	
	cout << "> diffusion step: " << steps << " sum: " << sum << " sum norm: "
	     << sum/(N*N) << " diff: " << diff << endl;
	

	OO(steps,0) = steps;
	OO(steps,1) = sum;
	OO(steps,2) = sum/(N*N);
	OO(steps,3) = diff;


	// Normalize each new obtained networks.
	cout << "run normaliseP: " << endl;
	cout << "P1new "; model->normaliseP( P1new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P2new "; model->normaliseP( P2new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "done." << endl; 
	
	cout << "run copyP... ";    
	model->copyP( P1new, P2new, model->PrDiag1, model->PrDiag2, model->map2Arr, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;    
	
    
	}//for
  
    

  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;
    runPcNorm = false;
  }
 

  //free memory
  model->freePrDiagSpace();
  model->freeSrSpace();
  
  if( runPcNorm ){
  
    try{
    
      // construct the combined affinity matrix by summing diffused matrices
      cout << "run updatePc... "; 
      model->updatePc( Pc, model->map2Arr, P1new, P2new, sum );
      Rcpp::checkUserInterrupt();
      cout << " done." << endl;
      
      // ensure affinity matrix is normalised and symmetrical
      ////cout << "run normalisePc: " << endl;
      ////cout << "Pc "; normalisePc( Pc, N, sum ); cout << sum << " ...done." << endl;
      
    } catch (Rcpp::internal::InterruptedException& e)
      {
	Rcout << "Caught an interrupt!" << std::endl;
      }
    
  }//if

  //free mapping vector
  model->freeMapSpace();
  
  
  //free Pnew space 
  if(P1new){free(P1new);}
  if(P2new){free(P2new);}
  
  //check we've freed space
  model->freeSpace();
   
  
  return Rcpp::List::create(Rcpp::Named("RESULT")=OO);

  
}



// [[Rcpp::export]]  
Rcpp::List fuseASDsets(Environment Pfuse, IntegerVector nCORES,
			 IntegerVector nTHREADS=1, IntegerVector nSTEPS=20 ){

    
  // pointer to Pc matrix object on disk
  XPtr<FBM_RW> pc = Pfuse["address_rw"];
  BMAcc_RW<double> Pc(pc);

  Rcpp::NumericMatrix OO;
  
  if(model==nullptr){
    cout << "create model first." << endl;
    OO = Rcpp::NumericMatrix(0,4);
    return Rcpp::List::create(Rcpp::Named("RESULT")=OO);
  }
  
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //--- Set nTHREADS
  if( nTHREADS[0] == 1 ){ nTHREADS[0] = nCORES[0]; }

  std::cout << "nTHREADS: " << nTHREADS[0] << endl;
  
  //---Setup fusion convergous 
  double FUSE_ACC     = 1.0e-6;
  int    FUSE_MAXSTEP = nSTEPS[0]+1;
  int    steps        = 0; 
  double delta        = 0;
  double delta_old    = 0;
  double diff         = 1;
  double sum          = 0;
  bool runPcNorm      = true;  

  std::cout << "Max number of steps: " << FUSE_MAXSTEP << endl;
  
  //set output
  OO  = Rcpp::NumericMatrix((FUSE_MAXSTEP+1),4);
  for(int i=0; i<(FUSE_MAXSTEP+1); i++){
    OO(i,0) = 0;
    OO(i,1) = 0;
    OO(i,2) = 0;
    OO(i,3) = 0;
  } 
  
  //set matrix size
  uint_fast32_t N = Pc.nrow(); 

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  if( model->map2Arr.empty() ){
  
    //mapping vector, make sure we reserve correct amount of space.
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K);

    //get a mapping from 2d matrix (with diagonal) to 1d array
    model->mapMatrixDiag2Array( N, model->map2Arr ); 

  }
  
  double *P1new = (double*)malloc(K*sizeof(double));
  double *P2new = (double*)malloc(K*sizeof(double));
  double *P3new = (double*)malloc(K*sizeof(double));
  double *P4new = (double*)malloc(K*sizeof(double));
 
  print_mem_usage();
  
  try{
     
      // check convergence 
      cout << "run updatePc... ";
      model->convergePc( model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, sum, N );
      Rcpp::checkUserInterrupt();
      cout << "done." << endl;     
      
      delta = sum;
      diff  = fabs(delta-delta_old);
  
      cout << "> diffusion step: " << 0 << " sum: " << sum << " sum norm: " << sum/(N*N)
	   << endl;

      OO(0,0) = 0;
      OO(0,1) = sum;
      OO(0,2) = sum/(N*N);
      OO(0,3) = diff;

    
      // perform the diffusion for t iterations
      cout << "perform diffusion:" << endl;
      for( steps=1; steps <= FUSE_MAXSTEP; steps++ ){
  
	// update each dataset and normalize each new obtained dataset.
	cout << "run updateP:" << endl;	
	cout << "P1new "; model->updateP( P1new, N, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->Sr1csc, model->tSr1csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P2new "; model->updateP( P2new, N, model->PrDiag1, model->PrDiag3, model->PrDiag4, model->Sr2csc, model->tSr2csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P3new "; model->updateP( P3new, N, model->PrDiag1, model->PrDiag2, model->PrDiag4, model->Sr3csc, model->tSr3csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "P4new "; model->updateP( P4new, N, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->Sr4csc, model->tSr4csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "done." << endl;

	// check convergence 
	cout << "run convergePc... ";
	model->convergePc( P1new, P2new, P3new, P4new, model->map2Arr, sum, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;     

	delta_old = delta;
	delta     = sum;
	diff      = fabs(delta-delta_old);
	
	cout << "> diffusion step: " << steps << " sum: " << sum << " sum norm: "
	     << sum/(N*N) << " diff: " << diff << endl;
	

	OO(steps,0) = steps;
	OO(steps,1) = sum;
	OO(steps,2) = sum/(N*N);
	OO(steps,3) = diff;

	cout << "run copyP... ";    
	model->copyP( P1new, P2new, P3new, P4new, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;    
	
    
	}//for
  
    

  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;
    runPcNorm = false;
  }
 

  //free memory 
  model->freePrDiagSpace();
  model->freeCSFSpace();
  
  if( runPcNorm ){
  
    try{
    
      // construct the combined affinity matrix by summing diffused matrices
      cout << "run updatePc... "; 
      model->updatePc( Pc, model->map2Arr, P1new, P2new, P3new, P4new, sum );
      Rcpp::checkUserInterrupt();
      cout << " done." << endl;     
       
    } catch (Rcpp::internal::InterruptedException& e)
      {
	Rcout << "Caught an interrupt!" << std::endl;
      }
    
  }//if

  cout << "flg0" << endl;
  
  //free mapping vector  
  model->freeMapSpace();  

  cout << "flg1" << endl;
  
  //free Pnew space 
  if(P1new){free(P1new);}
  if(P2new){free(P2new);}
  if(P3new){free(P3new);}
  if(P4new){free(P4new);}

  cout << "flg2" << endl;
  
  //check we've freed space
  //model->freeSpace();
 

  return Rcpp::List::create(Rcpp::Named("RESULT")=OO);

  
}


// [[Rcpp::export]]  
Rcpp::List fuseCONTROLsets(Environment Pfuse, IntegerVector nCORES,
			   IntegerVector nTHREADS=1, IntegerVector nSTEPS=20 ){

    
  // pointer to Pc matrix object on disk
  XPtr<FBM_RW> pc = Pfuse["address_rw"];
  BMAcc_RW<double> Pc(pc);

  Rcpp::NumericMatrix OO;
  
  if(model==nullptr){
    cout << "create model first." << endl;
    OO = Rcpp::NumericMatrix(0,4);
    return Rcpp::List::create(Rcpp::Named("RESULT")=OO);
  }
  
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] );

  //--- Set nTHREADS
  if( nTHREADS[0] == 1 ){ nTHREADS[0] = nCORES[0]; }

  std::cout << "nTHREADS: " << nTHREADS[0] << endl;
  
  //---Setup fusion convergous 
  double FUSE_ACC     = 1.0e-6;
  int    FUSE_MAXSTEP = nSTEPS[0]+1;
  int    steps        = 0; 
  double delta        = 0;
  double delta_old    = 0;
  double diff         = 1;
  double sum          = 0;
  bool runPcNorm      = true;  

  std::cout << "Max number of steps: " << FUSE_MAXSTEP << endl;
  
  //set output
  OO  = Rcpp::NumericMatrix((FUSE_MAXSTEP+1),4);
  for(int i=0; i<(FUSE_MAXSTEP+1); i++){
    OO(i,0) = 0;
    OO(i,1) = 0;
    OO(i,2) = 0;
    OO(i,3) = 0;
  } 
  
  //set matrix size
  uint_fast32_t N = Pc.nrow(); 

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  if( model->map2Arr.empty() ){
  
    //mapping vector, make sure we reserve correct amount of space.
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K);

    //get a mapping from 2d matrix (with diagonal) to 1d array
    model->mapMatrixDiag2Array( N, model->map2Arr ); 

  }
  
  double *P1new = (double*)malloc(K*sizeof(double));
  double *P2new = (double*)malloc(K*sizeof(double)); 
 
  print_mem_usage();
  
  try{
     
      // check convergence 
      cout << "run updatePc... ";
      model->convergePc( model->PrDiag1, model->PrDiag2, model->map2Arr, sum, N );
      Rcpp::checkUserInterrupt();
      cout << "done." << endl;     
      
      delta = sum;
      diff  = fabs(delta-delta_old);
  
      cout << "> diffusion step: " << 0 << " sum: " << sum << " sum norm: " << sum/(N*N)
	   << endl;

      OO(0,0) = 0;
      OO(0,1) = sum;
      OO(0,2) = sum/(N*N);
      OO(0,3) = diff;

    
      // perform the diffusion for t iterations
      cout << "perform diffusion:" << endl;
      for( steps=1; steps <= FUSE_MAXSTEP; steps++ ){
  
	// update each dataset and normalize each new obtained dataset.
	cout << "run updateP:" << endl;	
	cout << "P1new "; model->updateP( P1new, N, model->PrDiag2, model->Sr1csc, model->tSr1csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P2new "; model->updateP( P2new, N, model->PrDiag1, model->Sr2csc, model->tSr2csr, model->map2Arr, nTHREADS[0] ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();       
	cout << "done." << endl;

	// check convergence 
	cout << "run convergePc... ";
	model->convergePc( P1new, P2new, model->map2Arr, sum, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;     

	delta_old = delta;
	delta     = sum;
	diff      = fabs(delta-delta_old);
	
	cout << "> diffusion step: " << steps << " sum: " << sum << " sum norm: "
	     << sum/(N*N) << " diff: " << diff << endl;
	

	OO(steps,0) = steps;
	OO(steps,1) = sum;
	OO(steps,2) = sum/(N*N);
	OO(steps,3) = diff;

	cout << "run copyP... ";    
	model->copyP( P1new, P2new, model->PrDiag1, model->PrDiag2, model->map2Arr, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;    
	
    
	}//for
  
    

  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;
    runPcNorm = false;
  }
 

  //free memory 
  model->freePrDiagSpace();
  model->freeCSFSpace();
  
  if( runPcNorm ){
  
    try{
    
      // construct the combined affinity matrix by summing diffused matrices
      cout << "run updatePc... "; 
      model->updatePc( Pc, model->map2Arr, P1new, P2new, sum );
      Rcpp::checkUserInterrupt();
      cout << " done." << endl;     
       
    } catch (Rcpp::internal::InterruptedException& e)
      {
	Rcout << "Caught an interrupt!" << std::endl;
      }
    
  }//if

  cout << "flg0" << endl;
  
  //free mapping vector  
  model->freeMapSpace();  

  cout << "flg1" << endl;
  
  //free Pnew space 
  if(P1new){free(P1new);}
  if(P2new){free(P2new);}

  cout << "flg2" << endl;
  
  //check we've freed space
  //model->freeSpace();
 

  return Rcpp::List::create(Rcpp::Named("RESULT")=OO);

  
}


/*
//// [[Rcpp::export]]  
Rcpp::List fuseDatasets4(Environment Pfuse,IntegerVector nCORES, IntegerVector nSTEPS=20 ){

    
  // pointer to Pc matrix object on disk
  XPtr<FBM_RW> pc = Pfuse["address_rw"];
  BMAcc_RW<double> Pc(pc);

  Rcpp::NumericMatrix OO;
  
  if(model==nullptr){
    cout << "create model first." << endl;
    OO = Rcpp::NumericMatrix(0,4);
    return Rcpp::List::create(Rcpp::Named("RESULT")=OO);
  }
  
  
  //--- Setup OpenMP
  model->setOpenMP( nCORES[0] ); 
  
  //---Setup fusion convergous 
  double FUSE_ACC     = 1.0e-6;
  int    FUSE_MAXSTEP = nSTEPS[0]+1;
  int    steps        = 0; 
  double delta        = 0;
  double delta_old    = 0;
  double diff         = 1;
  double sum          = 0;
  bool runPcNorm      = true;  

  std::cout << "Max number of steps: " << FUSE_MAXSTEP << endl;
  
  //set output
  OO = Rcpp::NumericMatrix((FUSE_MAXSTEP+1),4);
  for(int i=0; i<(FUSE_MAXSTEP+1); i++){
    OO(i,0) = 0;
    OO(i,1) = 0;
    OO(i,2) = 0;
    OO(i,3) = 0;
  } 
  
  //set matrix size
  uint_fast32_t N = Pc.nrow(); 

  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  if( model->map2Arr.empty() ){
  
    //mapping vector, make sure we reserve correct amount of space.
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K);

    //get a mapping from 2d matrix (with diagonal) to 1d array
    model->mapMatrixDiag2Array( N, model->map2Arr ); 

  }
  
  double *P1new = (double*)malloc(K*sizeof(double));
  double *P2new = (double*)malloc(K*sizeof(double));
  double *P3new = (double*)malloc(K*sizeof(double));
  double *P4new = (double*)malloc(K*sizeof(double));
 

  try{
     
      // check convergence 
      cout << "run updatePc... ";
      model->convergePc( model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, sum, N );
      Rcpp::checkUserInterrupt();
      cout << "done." << endl;     
      
      delta = sum;
      diff  = fabs(delta-delta_old);
  
      cout << "> diffusion step: " << 0 << " sum: " << sum << " sum norm: " << sum/(N*N)
	   << endl;

      OO(0,0) = 0;
      OO(0,1) = sum;
      OO(0,2) = sum/(N*N);
      OO(0,3) = diff;

    
      // perform the diffusion for t iterations
      cout << "perform diffusion:" << endl;
      for( steps=1; steps <= FUSE_MAXSTEP; steps++ ){
  
	// update each dataset and normalize each new obtained dataset.
	cout << "run updateP:" << endl;	
	cout << "P1new "; model->updateP( P1new, N, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->COOSr1, model->map2Arr ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P2new "; model->updateP( P2new, N, model->PrDiag1, model->PrDiag3, model->PrDiag4, model->COOSr2, model->map2Arr ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();    
	cout << "P3new "; model->updateP( P3new, N, model->PrDiag1, model->PrDiag2, model->PrDiag4, model->COOSr3, model->map2Arr ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "P4new "; model->updateP( P4new, N, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->COOSr4, model->map2Arr ); cout << " ...done." << endl;
	Rcpp::checkUserInterrupt();   
	cout << "done." << endl;

	// check convergence 
	cout << "run convergePc... ";
	model->convergePc( P1new, P2new, P3new, P4new, model->map2Arr, sum, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;     

	delta_old = delta;
	delta     = sum;
	diff      = fabs(delta-delta_old);
	
	cout << "> diffusion step: " << steps << " sum: " << sum << " sum norm: "
	     << sum/(N*N) << " diff: " << diff << endl;
	

	OO(steps,0) = steps;
	OO(steps,1) = sum;
	OO(steps,2) = sum/(N*N);
	OO(steps,3) = diff;

	// Normalize each new obtained networks.
	cout << "run normaliseP: " << endl;
	cout << "P1new "; model->normaliseP( P1new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P2new "; model->normaliseP( P2new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P3new "; model->normaliseP( P3new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "P4new "; model->normaliseP( P4new, model->map2Arr, N, sum ); cout << sum << " ...done." << endl;
	Rcpp::checkUserInterrupt();
	cout << "done." << endl; 
	
	
	cout << "run copyP... ";    
	model->copyP( P1new, P2new, P3new, P4new, model->PrDiag1, model->PrDiag2, model->PrDiag3, model->PrDiag4, model->map2Arr, N );
	Rcpp::checkUserInterrupt();
	cout << "done." << endl;    
	
    
	}//for
  
    

  } catch (Rcpp::internal::InterruptedException& e) {
    Rcout << "Caught an interrupt!" << std::endl;
    runPcNorm = false;
  }
 

  //free memory 
  model->freePrDiagSpace();
  model->freeSrCOOSpace();
  
  if( runPcNorm ){
  
    try{
    
      // construct the combined affinity matrix by summing diffused matrices
      cout << "run updatePc... "; 
      model->updatePc( Pc, model->map2Arr, P1new, P2new, P3new, P4new, sum );
      Rcpp::checkUserInterrupt();
      cout << " done." << endl;     
       
    } catch (Rcpp::internal::InterruptedException& e)
      {
	Rcout << "Caught an interrupt!" << std::endl;
      }
    
  }//if

  //free mapping vector  
  model->freeMapSpace();  
  
  //free Pnew space 
  if(P1new){free(P1new);}
  if(P2new){free(P2new);}
  if(P3new){free(P3new);}
  if(P4new){free(P4new);}

  
  //check we've freed space
  model->freeSpace();
 

  return Rcpp::List::create(Rcpp::Named("RESULT")=OO);

  
}
*/
