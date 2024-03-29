/**
 * @file	lcmatrix.cxx
 *
 * @brief	Unit tests for Lucee::Matrix class
 */

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
// eigen inlcudes
#include <Eigen/Core>

// std includes
#include <cmath>
#include <ctime>
#include <vector>

// This weirdness is needed as CLAPACK defines integer as long int but
// as far I can see LAPACK itself only uses int. This can be a
// potential cause for problems but I do not know how to fix in
// general. (Ammar Hakim Wed Feb 1 2012).
#ifdef HAVE_CLAPACKCMAKE
#include <f2c.h>
#include <clapack.h>
#else
typedef int integer;
#include <LcLapackDeclarations.h>
#endif

/**
 * Compute matrix-vector multiply. Output vector must be
 * pre-allocated. Note that the computation performed is
 *
 * out = m*mat*vec + v*out
 *
 * @param m Factor for accumulation.
 * @param mat Matrix for multiplication.
 * @param vec Vector for multiplication.
 * @param v Factor for accumulation.
 * @param out On output, holds the product.
 */
void 
matVec(double m, const Lucee::Matrix<double>& mat, const double* vec, double v, double *out)
{
  double tv;
  unsigned rows = mat.numRows(), cols = mat.numColumns();
  for (unsigned i=0; i<rows; ++i)
  {
    tv = 0.0;
    for (unsigned j=0; j<cols; ++j)
      tv += mat(i,j)*vec[j];
    out[i] = m*tv + v*out[i];
  }
}

void 
matVec2(double m, unsigned rows, unsigned cols,
  const std::vector<double>& mat, const double* vec, double v, double *out)
{
  double tv;
  for (unsigned i=0; i<rows; ++i)
  {
    tv = 0.0;
    for (unsigned j=0; j<cols; ++j)
      tv += mat[i*rows+j]*vec[j];
    out[i] = m*tv + v*out[i];
  }
}

void 
matVec3(double m, unsigned rows, unsigned cols,
  const std::vector<double>& mat, const double* vec, double v, double *out)
{
  double tv;
  const double *mm = &mat[0];
  for (unsigned i=0; i<rows; ++i)
  {
    const double *vv = &vec[0];
    tv = 0.0;
    for (unsigned j=0; j<cols; ++j)
      tv += (*mm++)*(*vv++);
    out[i] = m*tv + v*out[i];
  }
}

// y = alpha*A*x + beta*y
// out = m*mat*vec + v*out
void
matVec4(double m, int M, int N, const double *mat, const double* vec, double v, double *out)
{
  int LDA = M;
  int INCX = 1;
  int INCY = 1;
  char TRANS[] = "N";

// call BLAS routine to do the multiplication  
  dgemv_(TRANS, &M, &N, &m, (double*) mat, &LDA, (double*)&vec[0], &INCX, &v, &out[0], &INCY);
}

int
doMatVec(unsigned N, int p, int nloop, unsigned mvtype)
{
  int lower[1] = {0};
  int upper[1]; upper[0] = N;
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::Field<1, double> elcFieldIn(rgn, p, 1.0);
  Lucee::Field<1, double> elcFieldOut(rgn, p, 1.0);
  Lucee::ConstFieldPtr<double> fInPtr = elcFieldIn.createConstPtr();
  Lucee::FieldPtr<double> fOutPtr = elcFieldOut.createPtr();

// allocate a matrix
  Lucee::Matrix<double> m(p,p);
  m = 1.0;

  std::vector<double> masv(p*p);
  for (unsigned i=0; i<p; ++i)
    for (unsigned j=0; j<p; ++j)
      masv[p*i+j] = 1.0;

  Eigen::MatrixXd matEig(p,p);
  for (unsigned i=0; i<p; ++i)
    for (unsigned j=0; j<p; ++j)
      matEig(i,j) = 1.0;
  Eigen::VectorXd vecInEig(p), vecOutEig(p);

  for (unsigned nl=0; nl<nloop; ++nl)
  {
    for (unsigned i=0; i<N; ++i)
    {
      elcFieldIn.setPtr(fInPtr, i);
      elcFieldOut.setPtr(fOutPtr, i);
// do matrix/vector multiply
      if (mvtype==1)
        matVec(1.0, m, &fInPtr[0], 0.0, &fOutPtr[0]);
      else if (mvtype==2)
        matVec2(1.0, p, p, masv, &fInPtr[0], 0.0, &fOutPtr[0]);
      else if (mvtype==3)
        matVec3(1.0, p, p, masv, &fInPtr[0], 0.0, &fOutPtr[0]);
      else if (mvtype==4)
        matVec4(1.0, p, p, &m.first(), &fInPtr[0], 0.0, &fOutPtr[0]);
      else if (mvtype==5)
      {
        vecOutEig.noalias() = matEig*vecInEig;
        for (unsigned i=0; i<p; ++i) fOutPtr[i] = vecOutEig[p];
      }
    }
  }

  unsigned nops = p*p; // number of operations
  return nops;
}

int
doMatVecAsMatMat(unsigned N, int p, int nloop)
{
  int lower[1] = {0};
  int upper[1]; upper[0] = N;
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::Field<1, double> elcFieldIn(rgn, p, 1.0);
  Lucee::Field<1, double> elcFieldOut(rgn, p, 1.0);
  Lucee::ConstFieldPtr<double> fInPtr = elcFieldIn.createConstPtr();
  Lucee::FieldPtr<double> fOutPtr = elcFieldOut.createPtr();

// allocate a matrix
  Lucee::Matrix<double> m(p,p);
  m = 1.0;

  std::vector<double> masv(p*p);
  for (unsigned i=0; i<p; ++i)
    for (unsigned j=0; j<p; ++j)
      masv[p*i+j] = 1.0;

  Eigen::MatrixXd matEig(p,p);
  for (unsigned i=0; i<p; ++i)
    for (unsigned j=0; j<p; ++j)
      matEig(i,j) = 1.0;
  unsigned chunkSize = p;

  int nchunks = N/chunkSize, nrem = N % chunkSize;
  std::cout << "nchunks " << nchunks << " nrem " << nrem << std::endl;
  Eigen::MatrixXd fchunk(p,chunkSize), frem(nrem,nrem);
  Eigen::MatrixXd fchunkOut(p,chunkSize), fremOut(nrem,nrem);

  unsigned nmul = 0;  
  for (unsigned nl=0; nl<nloop; ++nl)
  {
    for (unsigned c=0; c<nchunks; ++c)
    {
      //std::cout << c*chunkSize << std::endl;      
      for (unsigned j=0; j<chunkSize; ++j)
      {
        for (unsigned i=0; i<p; ++i)
          fchunk(i,j) = elcFieldIn(c*chunkSize+j,i);
      }
// multiply to compute vol integral
      fchunkOut = matEig*fchunk;
      nmul = nmul+1;

      for (unsigned j=0; j<chunkSize; ++j)
        for (unsigned i=0; i<p; ++i)
          elcFieldOut(c*chunkSize+j,i) = fchunkOut(i,j);
    }
  }
  std::cout << "Num of mat-mat muls " << nmul/nloop << std::endl;  

  unsigned nops = p*p; // number of operations
  return nops;
}

int
checkImpl()
{
  Lucee::Matrix<double> m(3,3);
  unsigned c=1;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
      m(i,j) = c++;
  
  std::vector<double> matasv(3*3);
  c=1;
  for (unsigned i=0; i<3; ++i)
    for (unsigned j=0; j<3; ++j)
      matasv[3*i+j] = c++;

  c=1;
  std::vector<double> vin(3);
  for (unsigned i=0; i<3; ++i) vin[i] = c++;

  std::vector<double> vo1(3), vo2(3), vo3(3), vo4(3);
  matVec(1.0, m, &vin[0], 0.0, &vo1[0]);
  matVec2(1.0, 3, 3, matasv, &vin[0], 0.0, &vo2[0]);
  matVec3(1.0, 3, 3, matasv, &vin[0], 0.0, &vo3[0]);
  matVec4(1.0, 3, 3, &m.first(), &vin[0], 0.0, &vo4[0]);

  unsigned nFail = 0;
  for (unsigned i=0; i<3; ++i)
    if ((vo1[i]!=vo2[i]) || (vo2[i]!=vo3[i]) || (vo3[i]!=vo4[i]))
      nFail++;
  return nFail;
}

int
main(int argc, char **argv) 
{
  using namespace std;
  unsigned nops;

  unsigned nFail = checkImpl();
  std::cout << "Number of failed tests (should be 0) " << nFail << std::endl;
  clock_t start, end;
  
  // start = clock();
  // nops = doMatVec(50, 4, 10000, 1);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=4 [Simple]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 8, 10000, 1);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=8 [Simple]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 27, 10000, 1);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=27 [Simple]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 100, 1000, 1);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=100 [Simple]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 200, 1000, 1);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=200 [Simple]"  << endl;

  // std::cout << "--------" << std::endl;
  // start = clock();
  // nops = doMatVec(50, 4, 10000, 3);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=4 [++Indexing]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 8, 10000, 3);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=8 [++Indexing]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 27, 10000, 3);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=27 [++Indexing]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 100, 1000, 3);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=100 [++Indexing]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 200, 1000, 3);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=200 [++Indexing]"  << endl;  

  // std::cout << "--------" << std::endl;
  // start = clock();
  // nops = doMatVec(50, 4, 10000, 4);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=4 [BLAS]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 8, 10000, 4);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=8 [BLAS]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 27, 10000, 4);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=27 [BLAS]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 100, 1000, 4);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=100 [BLAS]"  << endl;

  // start = clock();
  // nops = doMatVec(50, 200, 1000, 4);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=200 [BLAS]"  << endl;

  start = clock();
  nops = doMatVec(1200, 100, 1000, 4);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=100 [BLAS]"  << endl;

  start = clock();
  nops = doMatVec(1200, 100, 1000, 5);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=100 [EIGEN]"  << endl;

  start = clock();
  nops = doMatVec(1200, 100, 1000, 1);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=100 [matVec]"  << endl;

  std::cout << "--------" << std::endl;
  // start = clock();
  // nops = doMatVecAsMatMat(50, 4, 10000);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=4 [BLAS: gemm]"  << endl;

  // start = clock();
  // nops = doMatVecAsMatMat(1200, 100, 1000);
  // end = clock();
  // cout << right << fixed << difftime(end,start)  << " p=100 [BLAS: gemm]"  << endl;
}
