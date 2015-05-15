/**
 * @file	lcmatrix.cxx
 *
 * @brief	Unit tests for Lucee::Matrix class
 */

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>

// std includes
#include <cmath>
#include <ctime>

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

int
doMatVec(unsigned N, int p, int nloop)
{
  int lower[2] = {0, 0};
  int upper[2];
  upper[0] = upper[1] = N;
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcFieldIn(rgn, p, 1.0);
  Lucee::Field<2, double> elcFieldOut(rgn, p, 1.0);
  Lucee::ConstFieldPtr<double> fInPtr = elcFieldIn.createConstPtr();
  Lucee::FieldPtr<double> fOutPtr = elcFieldOut.createPtr();

// allocate a matrix
  Lucee::Matrix<double> m(p,p);
  m = 1.0;

  for (unsigned nl=0; nl<nloop; ++nl)
  {
    for (unsigned i=0; i<N; ++i)
      for (unsigned j=0; j<N; ++j)
      {
        elcFieldIn.setPtr(fInPtr, i, j);
        elcFieldOut.setPtr(fOutPtr, i, j);
// do matrix/vector multiply
        matVec(1.0, m, &fInPtr[0], 0.0, &fOutPtr[0]);
      }
  }

  return 0;
}

int
main(int argc, char **argv) 
{
  using namespace std;

  clock_t start = clock();
  doMatVec(100, 4, 10000);
  clock_t end = clock();
  cout << right << fixed << difftime(end,start)  << " (100x100 p=4)" << endl;

  start = clock();
  doMatVec(100, 8, 10000);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (100x100 p=9)" << endl;

  start = clock();
  doMatVec(50, 9, 10000);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (50x50 p=9)" << endl;

  start = clock();
  doMatVec(25, 16, 10000);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " (25x25 p=16)" << endl;

  return 0;
}
