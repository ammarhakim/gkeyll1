/**
 * @file	lceigen.cxx
 *
 * @brief	Unit tests for using EIGEN in Lucee
 */

// lucee includes
#include <LcLinAlgebra.h>
#include <LcTest.h>

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
    { return data[i+j*r]; }

    inline
    double& operator()(unsigned i, unsigned j)
    { return data[i+j*r]; }

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
      c(i,j) = 0.0;
    }

  clock_t start_t = clock();
// compute c
  for (unsigned i=0; i<N; ++i)
    for (unsigned j=0; j<N; ++j)
      for (unsigned k=0; k<N; ++k)
        c(i,j) += a(i,k)*b(k,j);
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

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lceigen");
  test_1();

  std::cout << "**** N = 500" << std::endl;
// N = 500
  double ht = handMM(500); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  double bt = blasMM(500); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  double et = eigenMM(500); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;

  std::cout << "\n**** N = 1000" << std::endl;
// N = 1000
  ht = handMM(1000); // HAND WRITTEN CODE
  std::cout << "Hand-written code took " << ht << "s" << std::endl;
  bt = blasMM(1000); // BLAS
  std::cout << "BLAS code took " << bt << "s" << std::endl;
  et = eigenMM(1000); // EIGEN
  std::cout << "Eigen code took " << et << "s" << std::endl;

  LC_END_TESTS;
}
