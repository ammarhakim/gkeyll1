/**
 * @file	lcmatrix.cxx
 *
 * @brief	Unit tests for Lucee::Matrix class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcMatrix.h>
#include <LcVector.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  Lucee::Matrix<double> A(2,5);

  LC_ASSERT("Testing if numRows is correct", A.numRows() == 2);
  LC_ASSERT("Testing if numColumns is correct", A.numColumns() == 5);

  unsigned shape[2] = {2,5};
  Lucee::Matrix<double> B(shape);

  LC_ASSERT("Testing if numRows is correct", B.numRows() == 2);
  LC_ASSERT("Testing if numColumns is correct", B.numColumns() == 5);

  Lucee::Matrix<double> S(2,2);
  S(0,0) = 1.0; S(0,1) = 0.0;
  S(1,0) = 0.0; S(1,1) = 2.0;

  Lucee::Vector<double> evr(2);
  Lucee::Vector<double> evi(2);

// compute eigenvalues of matrix
  S.eig(evr, evi);
  LC_ASSERT("Testing eigenvalues of S", evr[0]==1.0);
  LC_ASSERT("Testing eigenvalues of S", evr[1]==2.0);

  S(0,0) = 1.0; S(0,1) = 8.0;
  S(1,0) = 8.0; S(1,1) = 1.0;

  S.eig(evr, evi);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

// compute eigenvalues and eigenvectors of matrix
  Lucee::Matrix<double> vecl(2,2), vecr(2,2);
  S.eig(evr, evi, vecl, vecr);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

// check A*r = lambda*r for all right eigenvectors
  Lucee::Vector<double> eigVec(2);
  for (unsigned p=0; p<2; ++p)
  {
    eigVec = 0.0;
    for (unsigned i=0; i<2; ++i)
    {
      double sum = 0.0;
      for (unsigned j=0; j<2; ++j)
        sum += S(i,j)*vecr(j,p);
      eigVec[i] = sum;
    }
    for (unsigned i=0; i<2; ++i)
      LC_ASSERT("Checking A*r=lambda*r", epsCmp(eigVec[i], evr[p]*vecr(i,p)));
  }

// compute eigenvalues and left-eigenvectors of matrix
  Lucee::Matrix<double> vec(2,2);
  S.eigLeft(evr, evi, vec);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

  S.eigRight(evr, evi, vec);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

// check A*r = lambda*r for all right eigenvectors
  for (unsigned p=0; p<2; ++p)
  {
    eigVec = 0.0;
    for (unsigned i=0; i<2; ++i)
    {
      double sum = 0.0;
      for (unsigned j=0; j<2; ++j)
        sum += S(i,j)*vecr(j,p);
      eigVec[i] = sum;
    }
    for (unsigned i=0; i<2; ++i)
      LC_ASSERT("Checking A*r=lambda*r", epsCmp(eigVec[i], evr[p]*vecr(i,p)));
  }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcmatrix");
  test_1();
  LC_END_TESTS;
}
