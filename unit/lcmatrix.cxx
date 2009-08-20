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
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcmatrix");
  test_1();
  LC_END_TESTS;
}
