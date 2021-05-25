#include "Headers.h"

#ifndef FUSE_H
#define FUSE_H



//*** STRUCTS ***

//struct to store max delta
struct RESULT { double param1; };

//to use on row-major index ordering
struct sortRowCOO{
  bool operator()(const std::tuple<uint_fast32_t, uint_fast32_t, double> &a,
		  const std::tuple<uint_fast32_t, uint_fast32_t, double> &b) const {

    //std::tie, compare a[0] < b[0] and then a[1] < b[1]
    return std::tie(std::get<0>(a), std::get<1>(b)) < std::tie(std::get<0>(b), std::get<1>(a));    
   
    }
};

//to use on column-major index ordering
struct sortColCOO{
  bool operator()(const std::tuple<uint_fast32_t, uint_fast32_t, double> &a,
		  const std::tuple<uint_fast32_t, uint_fast32_t, double> &b) const {

    //std::tie, compare b[1] > a[1] and then b[0] > a[0]
    return std::tie(std::get<1>(b), std::get<0>(b)) > std::tie(std::get<1>(a), std::get<0>(a));
    
    }
};

//store coordinates, and weight, for an edge
typedef std::tuple<double, double, double> EDGE;

//store mapping between 2d matrix and 1d array
typedef std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t> tripleInt;

//store non-zero values in Sr matrcies, and index in 1d array - mapping a 2d diagonal matrix.
typedef std::tuple<uint_fast32_t, double, double> tripleSr;

//store non-zero values in Sr matrcies, and index in 1d array - mapping a 2d diagonal matrix.
typedef std::tuple<uint_fast32_t, uint_fast32_t, uint_fast32_t, double, double> tupleSr;

//Coordinate list (COO), store a list of (row, column, value) tuples.
typedef std::tuple<uint_fast32_t, uint_fast32_t, double> tupleCOO;


typedef struct CSF{

  arma::uword n_elem;
  arma::uword ncols;
  arma::uword nrows;
  arma::uword nz;
  
  double      *val;
  arma::uword *col;
  arma::uword *row;

  CSF() :
    n_elem( 0 ),
    ncols ( 0 ),
    nrows ( 0 ),
    nz    ( 0 ),
    val   ( nullptr ),
    col   ( nullptr ),
    row   ( nullptr ){}

  ~CSF(){
    freeSpace();
  }
 
  void freeSpace(){
       
    if(val){ delete[] val; }
    if(col){ delete[] col; }
    if(row){ delete[] row; }
    
    n_elem=0;
    ncols =0;
    nrows =0;
    nz    =0;
   
  }

  void init_csc( arma::uword NVALS, arma::uword NCOLS ){

    freeSpace();    

    arma::uword NCOLSP1 = static_cast<arma::uword>(NCOLS+1);
    
    //initialise CSC struct   
    val    = new double[NVALS];
    row    = new arma::uword[NVALS];
    col    = new arma::uword[NCOLSP1];
    nz     = NVALS;
    n_elem = NVALS;
    nrows  = NVALS;
    ncols  = NCOLS;
   

  }

  void zero_csc(){

    arma::uword i;
    
    if(val){ for(i=0; i<nz; i++)val[i] = 0; }
    if(row){ for(i=0; i<nz; i++)row[i] = 0; }
    if(col){ for(i=0; i<(ncols+1); i++)col[i] = 0; }
    
  }
  
  void init_csr( arma::uword NVALS, arma::uword NROWS ){

    freeSpace(); 
   
    arma::uword NROWSP1 = static_cast<arma::uword>(NROWS+1);
    
    //initialise CSR struct   
    val    = new double[NVALS];
    col    = new arma::uword[NVALS];
    row    = new arma::uword[NROWSP1];
    nz     = NVALS;
    n_elem = NVALS;
    ncols  = NVALS;
    nrows  = NROWS;
   
    
  }

  void zero_csr(){

    arma::uword i;
    
    if(val){ for(i=0; i<nz; i++)val[i] = 0; }
    if(col){ for(i=0; i<nz; i++)col[i] = 0; }
    if(row){ for(i=0; i<(nrows+1); i++)row[i] = 0; }
    
  }
  
} CSF;


//***************


//*** INLINES ***
inline void          int2fast32(int i, uint_fast32_t &i32){ i32 = static_cast<uint_fast32_t>(i); }
inline uint_fast32_t int2fast32(int i)                    { return static_cast<uint_fast32_t>(i); }

inline void fast32toint(uint_fast32_t i32, int &i){ i = static_cast<int>(i32); }
inline int  fast32toint(uint_fast32_t i32)        { return static_cast<int>(i32); }

inline float  double2float( double d ){ return static_cast<float>(d);  }
inline double float2double( float f  ){ return static_cast<double>(f); }


//Indexing from 2d matrix (with diagonal) to 1d array
inline uint_fast32_t Trag_eq(const uint_fast32_t row, const uint_fast32_t col, const uint_fast32_t N){
  //Reference:
  //1) https://www.codeguru.com/cpp/cpp/algorithms/general/article.php/c11211/TIP-Half-Size-Triangular-Matrix.htm
  //2) http:www.research.att.com/%7Enjas/sequences/A000217
  if (row <= col)
    return row*N - (row-1)*((row-1) + 1)/2 + col - row;
  else if (col<row)
    return col*N - (col-1)*((col-1) + 1)/2 + row - col;

  return -1;
}

//Indexing from 2d matrix (with-out diagonal) to 1d array
inline uint_fast32_t Trag_noeq(const uint_fast32_t row, const uint_fast32_t col, const uint_fast32_t N){
  //assert(row != col);    //You can add this in if you like
  if (row<col)
    return row*(N-1) - (row-1)*((row-1) + 1)/2 + col - row - 1;
  else if (col<row)
    return col*(N-1) - (col-1)*((col-1) + 1)/2 + row - col - 1;
  else
    return -1;
}

//Indexing from 1d array to 2d matrix (with diagonal)
inline void Trag_reverse_eq(const uint_fast32_t index, const uint_fast32_t N, uint_fast32_t& row, uint_fast32_t& col) { 

  row = 0; 
  uint_fast32_t keyafter; 
  do 
    { 
      row++; 
      keyafter = row * N - (row - 1) * row / 2; 
    } while (index >= keyafter); 
  row--; 
  col = N - keyafter + index; 

} 

//Indexing from 1d array to 2d matrix (with-out diagonal)
inline void Trag_reverse_noeq(const uint_fast32_t index, const uint_fast32_t N, uint_fast32_t& row, uint_fast32_t& col){

  row = 0;
  uint_fast32_t keyafter;
  do
    {
      row++;
      keyafter = row * (N-1) - (row - 1) * row / 2;
    } while (index >= keyafter);
  row--;
  col = N - keyafter + index;

}


inline arma::mat FBM2arma( Rcpp::Environment BM ) {
  //https://github.com/privefl/bigstatsr/blob/master/src/arma-prod.cpp
  
  Rcpp::XPtr<FBM> xpBM = BM["address"];
  myassert(xpBM->matrix_type() == 8,
           "Mapping to arma::mat is available for 'double' FBMs only.");

  return arma::mat((double*)xpBM->matrix(), xpBM->nrow(), xpBM->ncol(), false);
  
}

inline arma::mat Arr2arma( double *x, const uint_fast32_t xlen, const uint_fast32_t N, const uint_fast32_t M ) {

  myassert( xlen == N*M, "array length not equal to nrow*ncol."); 

  return arma::mat(x, N, M, false);

}

inline void Arr2ArmaD( arma::mat &Mat, const double *x, const uint_fast32_t N, const vector<tripleInt> &map ){

  myassert( Mat.n_elem == N*N, "Matrix size not equal to N*N." );
  
  uint_fast32_t i,j,k,kk,K;

  K = map.size();

#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,kk) shared(x,Mat,map)
  for(kk=0;kk<K;kk++){
    i = std::get<0>(map[kk]);
    j = std::get<1>(map[kk]);
    k = std::get<2>(map[kk]);			     

    Mat(i,j) = x[kk];
    Mat(j,i) = x[kk];
    
  }   
  
}

inline void SrVec2ArmaD( arma::mat &Mat, const vector<tupleCOO> &map, const uint_fast32_t N ){

  myassert( Mat.n_elem == N*N, "Matrix size not equal to N*N." );
  
  uint_fast32_t i,j,k,kk,K;

  double Val;

  K = map.size();

#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,kk,Val) shared(Mat,map)
  for(kk=0;kk<K;kk++){
    i   = std::get<0>(map[kk]);
    j   = std::get<1>(map[kk]);
    Val = std::get<2>(map[kk]);			     

    Mat(i,j) = Val;
        
  }   
  
}


inline void SrVec2ArmaD( arma::SpMat<double> &Mat, const vector<tupleCOO> &map, const uint_fast32_t N ){

  myassert( Mat.n_elem == N*N, "Matrix size not equal to N*N." );
  
  uint_fast32_t i,j,k,kk,K;

  double Val;

  K = map.size();

#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,kk,Val) shared(Mat,map)
  for(kk=0;kk<K;kk++){
    i   = std::get<0>(map[kk]);
    j   = std::get<1>(map[kk]);
    Val = std::get<2>(map[kk]);			     

    Mat(i,j) = Val;    
    
  }   
  
}


inline void tSrVec( const vector<tupleCOO> &map, vector<tupleCOO> &tmap, const uint_fast32_t N ){

  uint_fast32_t i,j,k,kk,K;

  double Val;

  K = map.size();

  tmap.reserve(K);
  
  #pragma omp declare reduction (merge : std::vector<tupleCOO> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  
#pragma omp parallel for schedule(static) default(none) firstprivate(K) private(i,j,kk,Val) shared(map) reduction(merge: tmap)
  for(kk=0;kk<K;kk++){
    i   = std::get<0>(map[kk]);
    j   = std::get<1>(map[kk]);
    Val = std::get<2>(map[kk]);			     

    tmap.push_back( tupleCOO(j,i,Val) );
    
  }   
  
}

/////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
//https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
inline void process_mem_usage(double& vm_usage, double& resident_set){

  using std::ios_base;
  using std::ifstream;
  using std::string;

  vm_usage     = 0.0;
  resident_set = 0.0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	      >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	      >> utime >> stime >> cutime >> cstime >> priority >> nice
	      >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

  stat_stream.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0;
  resident_set = rss * page_size_kb;

}

inline void print_mem_usage(){

  double vm, rss;
  process_mem_usage(vm, rss);
  std::cout << "*** MEMORY USEAGE ***" << std::endl;
  std::cout << "VM: " << vm << "; RSS: " << rss << std::endl;
  std::cout << "*********************" << std::endl;
}

//**************


class fuse {

protected: 

  
public:
   fuse();
  ~fuse();
  
  void setOpenMP( int, bool=true );
  
  void getOpenMPStart( double &);

  void getOpenMPEnd( double &);

  void freeCSFSpace();

  void freePrDiagSpace();

  void freePrBlckSpace();

  void freeSrSpace();

  void freeSrCOOSpace();

  void freeMapSpace();

  void setMapSpace(int);

  void freeSpace();

  void csr_dense_tcrossprod( const CSF *, const arma::Mat<double> &, arma::Mat<double> &, const int=1 );

  void dense_csc_prod(const arma::Mat<double>&, const CSF *, arma::Mat<double> &, const int=1 );

  void mappingCOO2CSF( vector<tupleCOO> &, CSF *, const uint_fast32_t,  bool=true );

  void mapSrmatrix(BMAcc_RW<double> &, const uint_fast32_t, vector<tupleCOO>&, const vector<tripleInt>&, bool=true);

  void mapSrmatrix(BMAcc_RW<double> &, const uint_fast32_t, vector<tripleSr>&, const vector<tripleInt>&);

  void mapSrmatrix( BMAcc_RW<double> &, const uint_fast32_t, CSF *, CSF *, const vector<tripleInt> & );

  void mapSrmatrix(BMAcc_RW<double> &, const uint_fast32_t, vector<tupleSr>&, const vector<tripleInt> &);

  void copyDatasetRtoC( double *, BMAcc_RW<double> &);

  void copyDatasetRtoC( BMAcc_RW<double> &, double *, const uint_fast32_t, const vector<tripleInt> &, double &, bool=true );

  void copyDatasetRtoC3( BMAcc_RW<double> &, double *, const uint_fast32_t, const vector<tripleInt> &, double &, bool=true );

  void copyDatasetCtoR( BMAcc_RW<double> &, const double * );

  void transpose( const double *, double *, const vector<tripleInt> &, const uint_fast32_t);

  void transpose( const double *, double *, const uint_fast32_t, bool);

  void transpose( const arma::Mat<double> &, arma::Mat<double> &, const uint_fast32_t, bool );

  void transpose( BMAcc_RW<double> &, double *, const uint_fast32_t, bool );

  void symmetrise( double *, const uint_fast32_t, double & );
  
  void symmetrise( BMAcc_RW<double> &, double *, const uint_fast32_t, double & );

  void symmetrise( arma::Mat<double> &, const uint_fast32_t, double & );

  void symmetrise( double *, const double *, const double *, const uint_fast32_t );

  void symmetrise( BMAcc_RW<double> &, const double *, const double *, const uint_fast32_t );

  void rowNorm( double *, const uint_fast32_t, bool=false );

  void initialiseP( BMAcc_RW<double> &, double *, const uint_fast32_t, const vector<tripleInt> &, double & );

  void updateP( double *, const uint_fast32_t, const double *, const double *, const double *, const vector<tupleSr> &, const vector<tripleInt> &, double & );

  void updateP( double *, const uint_fast32_t, const double *, const vector<tupleSr> & );
  
  void updateP( double *, const uint_fast32_t, const double *, const double *, const double *, const vector<tupleSr> & );

void updateP( double *, const uint_fast32_t, const double *, CSF *, CSF *, const vector<tripleInt> &, int=1 );
  
  void updateP( double *, const uint_fast32_t, const double *, const double *, const double *, CSF *, CSF *, const vector<tripleInt> &, int=1 );

  void updateP( double *, const uint_fast32_t, const double *, const double *, const double *, const vector<tupleCOO> &, const vector<tripleInt> & );

  void convergePc( const double *, const double *, const vector<tripleInt> &, double &, const uint_fast32_t );
  
  void convergePc( const double *, const double *, const double *, const double *, const vector<tripleInt> &, double &, const uint_fast32_t );


  void normaliseP( double *, const vector<tripleInt> &, const uint_fast32_t, double & );

  void normaliseP( double *, const uint_fast32_t, double & );

   void copyP( const double *, const double *, double *, double *,
	       const vector<tripleInt> &, const uint_fast32_t );
  
  void copyP( const double *, const double *, const double *, const double *,
	      double *, double *, double *, double *, const vector<tripleInt> &, const uint_fast32_t );

  void copyP( const double *, const double *, double *, double *, const uint_fast32_t );
  
  void copyP( const double *, const double *, const double *, const double *, double *, double *, double *, double *, const uint_fast32_t );

  void updatePc( BMAcc_RW<double> &, const vector<tripleInt> &, const double *, const double *, double & );
  
  void updatePc( BMAcc_RW<double> &, const vector<tripleInt> &, const double *, const double *, const double *, const double *, double & );

  void normalisePc( BMAcc_RW<double> &, const uint_fast32_t, double & );
  
  void map2Array( const double *, double *, const uint_fast32_t, const vector<tripleInt> & );

  void mapMatrixDiag2Array( const uint_fast32_t, vector<tripleInt> & );
  
  void mapArray2MatrixDiag( const uint_fast32_t, vector<tripleInt> & );

  //variables
  double EPSILON;
  double logEPSILON;

  double OMP_start;
  double OMP_end;

  CSF *tSr1csr;
  CSF *Sr1csc;
  CSF *tSr2csr;
  CSF *Sr2csc;
  CSF *tSr3csr;
  CSF *Sr3csc;
  CSF *tSr4csr;
  CSF *Sr4csc;
  
  double *PrDiag1;
  double *PrDiag2;
  double *PrDiag3;
  double *PrDiag4;

  double *PrBlck1;
  double *PrBlck2;
  double *PrBlck3;
  double *PrBlck4;

  bool freedCSFSpace;
  bool freedPrDiagSpace;
  bool freedPrBlckSpace;
  
  vector<tupleCOO>  COOSr1;
  vector<tupleCOO>  COOSr2;
  vector<tupleCOO>  COOSr3;
  vector<tupleCOO>  COOSr4; 
  
  vector<tupleSr>  mapSr1;
  vector<tupleSr>  mapSr2;
  vector<tupleSr>  mapSr3;
  vector<tupleSr>  mapSr4; 
  
  
  vector<tripleInt> map2Arr;


};

#endif

