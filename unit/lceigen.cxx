/**
 * @file	lceigen.cxx
 *
 * @brief	Unit tests for using EIGEN in Lucee
 */

// lucee includes
#include <LcLinAlgebra.h>
#include <LcTest.h>

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

// eigen inlcudes
#include <Eigen/Core>

// std includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <vector>

void
test_1()
{
  Eigen::MatrixXd m(4,4);
  for (unsigned i=0; i<4; ++i)
    for (unsigned j=0; j<4; ++j)
      m(i,j) = 1.0;

  LC_ASSERT("Testing rows", m.rows() == 4);
  LC_ASSERT("Testing cols", m.cols() == 4);

  m.resize(3,5);
  LC_ASSERT("Testing rows", m.rows() == 3);
  LC_ASSERT("Testing cols", m.cols() == 5);
}

class Md
{
  public:
    Md(unsigned r, unsigned c)
      : r(r), c(c), data(r*c)
    {}

    inline
    double operator()(unsigned i, unsigned j) const
    { return data[i*c+j]; }

    inline
    double& operator()(unsigned i, unsigned j)
    { return data[i*c+j]; }

  private:
    unsigned r, c;
    std::vector<double> data;
};

double
handMM(unsigned N)
{
  Md a(N,N), b(N,N), c(N,N);

// fill a and b with random numbers
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
    {
      a(i,j) = std::rand()/RAND_MAX;
      b(i,j) = std::rand()/RAND_MAX;
    }

  clock_t start_t = clock();
// compute c
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
    {
      c(i,j) = 0.0;
      for (unsigned k=0; k<N; ++k)
        c(i,j) += a(i,k)*b(k,j);
    }
  clock_t end_t = clock();

  return (double) (end_t-start_t)/CLOCKS_PER_SEC;
}

double
handMM2(unsigned N)
{
  Md a(N,N), b(N,N), c(N,N), tmp(N,N);

// fill a and b with random numbers
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
    {
      a(i,j) = std::rand()/RAND_MAX;
      b(i,j) = std::rand()/RAND_MAX;
    }

  clock_t start_t = clock();

  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
      tmp(i,j) = b(j,i);

// compute c
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
    {
      c(i,j) = 0.0;
      for (unsigned k=0; k<N; ++k)
        c(i,j) += a(i,k)*tmp(j,k);
    }
  clock_t end_t = clock();

  return (double) (end_t-start_t)/CLOCKS_PER_SEC;
}

double
eigenMM(unsigned N)
{
  Eigen::MatrixXd a(N,N), b(N,N), c(N,N);

// fill a and b with random numbers
  for (unsigned j=0; j<N; ++j)
    for (unsigned i=0; i<N; ++i)
    {
      a(i,j) = std::rand()/RAND_MAX;
      b(i,j) = std::rand()/RAND_MAX;
    }
// compute c
  clock_t start_t = clock();
  c = a*b;
  clock_t end_t = clock();

  return (double) (end_t-start_t)/CLOCKS_PER_SEC;
}

double
blasMM(unsigned N)
{
  Lucee::Matrix<double> a(N,N), b(N,N), c(N,N);

// fill a and b with random numbers
  for (unsigned j=0; j<N; ++j)
    for (unsigned i=0; i<N; ++i)
    {
      a(i,j) = std::rand()/RAND_MAX;
      b(i,j) = std::rand()/RAND_MAX;
    }
// compute c
  clock_t start_t = clock();
  Lucee::accumulate(0.0, c, 1.0, a, b);
  clock_t end_t = clock();

  return (double) (end_t-start_t)/CLOCKS_PER_SEC;
}

// c = a*b
void
blasMatMat(int shapea[2], double *a, int shapeb[2], double *b, int shapec[2], double *c)
{
  char TRANSA[] = "N", TRANSB[] = "N";
  int M = shapea[0], N = shapeb[1];
  int K = shapeb[1]; // shapeb[0];
  double alpha = 1.0, beta = 0.0;
  int LDA = M, LDB = K, LDC = M;
  dgemm_(TRANSA, TRANSB, &M, &N, &K, &alpha, a, &LDA, b, &LDB, &beta, c, &LDC);
}

void
blasMatVec(int shapea[2], double *a, double *x, double *y)
{
  char TRANS[] = "N";
  int M = shapea[0], N = shapea[1];
  double alpha = 1.0, beta = 0.0;
  int LDA = M, INCX=1, INCY=1;
  dgemv_(TRANS, &M, &N, &alpha, a, &LDA, x, &INCX, &beta, y, &INCY);
}

Eigen::MatrixXd
eigenMatMat(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
{
  return a*b;
}

int
test_mat_mat_vs_mat_vec(unsigned N, unsigned ntries)
{
  int shape[2];
  shape[0] = N; shape[1] = N;
  unsigned nelem = shape[0]*shape[1];
  std::vector<double> a(nelem), b(nelem), c(nelem);
  std::vector<double> x(N), y(N);  

// initialize stuff
  for (unsigned i=0; i<nelem; ++i)
  {
    a[i] = std::rand()/RAND_MAX;
    b[i] = std::rand()/RAND_MAX;
  }
  for (unsigned i=0; i<N; ++i)
    x[i] = std::rand()/RAND_MAX;

  clock_t start_t, end_t;

  start_t = clock();
// time mat-mat  
  for (unsigned n=0; n<ntries; ++n)
    blasMatMat(shape, &a[0], shape, &b[0], shape, &c[0]);
  end_t = clock();
  double tMM = (double) (end_t-start_t)/CLOCKS_PER_SEC;

  start_t = clock();  
// time mat-vec  
  for (unsigned n=0; n<ntries; ++n)
  {
    for (unsigned k=0; k<N; ++k)
      blasMatVec(shape, &a[0], &x[0], &y[0]);
  }
  end_t = clock();
  double tMV = (double) (end_t-start_t)/CLOCKS_PER_SEC;

  // std::cout << "Timing for (N, ntries) = (" << N << ", " << ntries << ") " << std::endl;
  // std::cout << "   Mat-Mat --> " << tMM << std::endl;
  // std::cout << "   Mat-Vec --> " << tMV << std::endl;
  // std::cout << std::endl;
  std::cout << N << " " << tMM << " " << tMV << std::endl;
}

int
test_mat_mat()
{
  double ht, htc, bt, et;

  std::cout << "**** N = 27" << std::endl;  
// N = 27
  ht = handMM(27); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  htc = handMM2(27); // HAND WRITTEN CODE
  std::cout << "Hand-written cache-friendly code took " << htc << "s" << std::endl;
  bt = blasMM(27); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  et = eigenMM(27); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;    

  std::cout << "**** N = 100" << std::endl;  
// N = 100
  ht = handMM(100); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  htc = handMM2(100); // HAND WRITTEN CODE
  std::cout << "Hand-written cache-friendly code took " << htc << "s" << std::endl;
  bt = blasMM(100); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  et = eigenMM(100); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;  

  std::cout << "**** N = 500" << std::endl;  
// N = 500
  ht = handMM(500); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  htc = handMM2(500); // HAND WRITTEN CODE
  std::cout << "Hand-written cache-friendly code took " << htc << "s" << std::endl;
  bt = blasMM(500); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  et = eigenMM(500); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;

  std::cout << "\n**** N = 1000" << std::endl;
// N = 1000
  ht = handMM(1000); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  htc = handMM2(1000); // HAND WRITTEN CODE
  std::cout << "Hand-written cache-friendly code took " << htc << "s" << std::endl;
  bt = blasMM(1000); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  et = eigenMM(1000); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;
}

int
main(int argc, char *argv[])
{
  //LC_BEGIN_TESTS("lceigen");
  //test_1();

  test_mat_mat_vs_mat_vec(25, 100);
  test_mat_mat_vs_mat_vec(50, 100);
  test_mat_mat_vs_mat_vec(75, 100);
  test_mat_mat_vs_mat_vec(100, 100);
  test_mat_mat_vs_mat_vec(150, 100);
  test_mat_mat_vs_mat_vec(200, 100);
  test_mat_mat_vs_mat_vec(250, 100);
  test_mat_mat_vs_mat_vec(300, 100);
  test_mat_mat_vs_mat_vec(400, 100);
  test_mat_mat_vs_mat_vec(500, 100);
  test_mat_mat_vs_mat_vec(750, 100);
  
  //LC_END_TESTS;
}
