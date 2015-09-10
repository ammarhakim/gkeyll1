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
#include <vector>

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

int
doMatVec(unsigned N, int p, int nloop, unsigned mvtype)
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

  std::vector<double> masv(p*p);
  for (unsigned i=0; i<p; ++i)
    for (unsigned j=0; j<p; ++j)
      masv[p*i+j] = 1.0;

  for (unsigned nl=0; nl<nloop; ++nl)
  {
    for (unsigned i=0; i<N; ++i)
      for (unsigned j=0; j<N; ++j)
      {
        elcFieldIn.setPtr(fInPtr, i, j);
        elcFieldOut.setPtr(fOutPtr, i, j);
// do matrix/vector multiply
        if (mvtype==1)
          matVec(1.0, m, &fInPtr[0], 0.0, &fOutPtr[0]);
        else if (mvtype==2)
          matVec2(1.0, p, p, masv, &fInPtr[0], 0.0, &fOutPtr[0]);
        else if (mvtype==3)
          matVec3(1.0, p, p, masv, &fInPtr[0], 0.0, &fOutPtr[0]);
      }
  }

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

  std::vector<double> vo1(3), vo2(3), vo3(3);
  matVec(1.0, m, &vin[0], 0.0, &vo1[0]);
  matVec2(1.0, 3, 3, matasv, &vin[0], 0.0, &vo2[0]);
  matVec3(1.0, 3, 3, matasv, &vin[0], 0.0, &vo3[0]);

  unsigned nFail = 0;
  for (unsigned i=0; i<3; ++i)
    if ((vo1[i]!=vo2[i]) || (vo2[i]!=vo3[i]))
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
  
  clock_t start = clock();
  nops = doMatVec(50, 4, 5000, 1);
  clock_t end = clock();
  cout << right << fixed << difftime(end,start)  << " p=4 [1]"  << endl;

  start = clock();
  nops = doMatVec(50, 8, 5000, 1);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=8 [1]"  << endl;

  start = clock();
  nops = doMatVec(50, 21, 5000, 1);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=21 [1]"  << endl;

  std::cout << "--------" << std::endl;
  start = clock();
  nops = doMatVec(50, 4, 5000, 2);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=4 [2]"  << endl;

  start = clock();
  nops = doMatVec(50, 8, 5000, 2);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=8 [2]"  << endl;

  start = clock();
  nops = doMatVec(50, 21, 5000, 2);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=21 [2]"  << endl;

  std::cout << "--------" << std::endl;
  start = clock();
  nops = doMatVec(50, 4, 5000, 3);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=4 [3]"  << endl;

  start = clock();
  nops = doMatVec(50, 8, 5000, 3);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=8 [3]"  << endl;

  start = clock();
  nops = doMatVec(50, 21, 5000, 3);
  end = clock();
  cout << right << fixed << difftime(end,start)  << " p=21 [3]"  << endl;
}
