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
#include <LcExcept.h>
#include <LcLinAlgebra.h>
#include <LcMatrix.h>
#include <LcTest.h>
#include <LcVector.h>

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
  Lucee::eig(S, evr, evi);
  LC_ASSERT("Testing eigenvalues of S", evr[0]==1.0);
  LC_ASSERT("Testing eigenvalues of S", evr[1]==2.0);

  S(0,0) = 1.0; S(0,1) = 8.0;
  S(1,0) = 8.0; S(1,1) = 1.0;

  Lucee::eig(S, evr, evi);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

// compute eigenvalues and eigenvectors of matrix
  Lucee::Matrix<double> vecl(2,2), vecr(2,2);
  Lucee::eig(S, evr, evi, vecl, vecr);
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

// check l*A = lambda*l for all left eigenvectors
  for (unsigned p=0; p<2; ++p)
  {
    eigVec = 0.0;
    for (unsigned j=0; j<2; ++j)
    {
      double sum = 0.0;
      for (unsigned i=0; i<2; ++i)
        sum += vecl(i,p)*S(i,j);
      eigVec[j] = sum;
    }
    for (unsigned i=0; i<2; ++i)
      LC_ASSERT("Checking l*A=l*lambda", epsCmp(eigVec[i], evr[p]*vecl(i,p)));
  }

// compute eigenvalues and left-eigenvectors of matrix
  Lucee::Matrix<double> vec(2,2);
  Lucee::eigLeft(S, evr, evi, vec);
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[0], 9.0));
  LC_ASSERT("Testing eigenvalues of S", epsCmp(evr[1], -7.0));

// check l*A = lambda*l for all left eigenvectors
  for (unsigned p=0; p<2; ++p)
  {
    eigVec = 0.0;
    for (unsigned j=0; j<2; ++j)
    {
      double sum = 0.0;
      for (unsigned i=0; i<2; ++i)
        sum += vecl(i,p)*S(i,j);
      eigVec[j] = sum;
    }
    for (unsigned i=0; i<2; ++i)
      LC_ASSERT("Checking l*A=l*lambda", epsCmp(eigVec[i], evr[p]*vecl(i,p)));
  }

// compute eigenvalues and right-eigenvectors of matrix
  Lucee::eigRight(S, evr, evi, vec);
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

void
test_2()
{
  unsigned shape[2] = {2, 5};
  int start[2] = {-1, -5};
  Lucee::Matrix<double> A(shape, start);
  A = 10.5;
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing matrix assignment operator", A(i,j) == 10.5);
}

void
test_3()
{
  unsigned shape[2] = {2, 5};
  int start[2] = {-1, -5};
  Lucee::Matrix<double> A(shape, start);
  Lucee::Matrix<double> ACopy(A);
  A = 10.5;
 
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing matrix assignment operator", ACopy(i,j) == 10.5);
  
  ACopy = 12.5;
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing matrix assignment operator", A(i,j) == 12.5);
}

void
test_4()
{
  unsigned shape[2] = {2, 5};
  int start[2] = {-1, -5};
  Lucee::Matrix<double> A(shape, start);

  unsigned ones[2] = {1, 1};
  int zeros[2] = {0, 0}; 
  Lucee::Matrix<double> AAss(ones, zeros);
  AAss = A;
 
  A = 10.5;
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing matrix assignment operator", AAss(i,j) == 10.5);
  
  AAss = 12.5;
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing matrix assignment operator", A(i,j) == 12.5);
}

void
test_5()
{
  Lucee::Matrix<double> S(2,3);
  Lucee::Vector<double> evr(2);
  Lucee::Vector<double> evi(2);
  Lucee::Matrix<double> vecr(2,2);
  Lucee::Matrix<double> vecl(2,2);

  LC_RAISES("Testing if exception is thrown", Lucee::eig(S, evr, evi), Lucee::Except);
  LC_RAISES("Testing if exception is thrown", Lucee::eig(S, evr, evi, vecl, vecr), Lucee::Except);
  LC_RAISES("Testing if exception is thrown", Lucee::eigRight(S, evr, evi, vecl), Lucee::Except);
  LC_RAISES("Testing if exception is thrown", Lucee::eigLeft (S, evr, evi, vecl), Lucee::Except);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcmatrix");
  test_1();
  test_2();
  test_3();
  test_4();
  test_5();
  LC_END_TESTS;
}
