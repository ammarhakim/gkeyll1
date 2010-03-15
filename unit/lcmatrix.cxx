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

void
test_6()
{
  Lucee::Matrix<double> S(2,2);
  S(0,0) = 1.0; S(0,1) = 2.0;
  S(1,0) = 3.0; S(1,1) = 4.0;

  Lucee::Matrix<double> Sdup(S.duplicate());
  LC_ASSERT("Testing duplicated matrix", Sdup(0,0) == 1.0);
  LC_ASSERT("Testing duplicated matrix", Sdup(0,1) == 2.0);
  LC_ASSERT("Testing duplicated matrix", Sdup(1,0) == 3.0);
  LC_ASSERT("Testing duplicated matrix", Sdup(1,1) == 4.0);
}

void
test_7()
{
  Lucee::Matrix<double> A(2,3), B(3,4), C(2,4);

  A = 1.0;
  B = 2.0;
  Lucee::accumulate(C, A, B); // C = A*B
// check if the accumulate worked
  for (int i=C.getLower(0); i<C.getUpper(0); ++i)
    for (int j=C.getLower(1); j<C.getUpper(1); ++j)
      LC_ASSERT("Testing if accumulate worked", C(i,j)==6.0);

  Lucee::Matrix<double> A1(2,2), B1(2,1), C1(2,1);
  A1(0,0) = 2.0; A1(0,1) = 1.0;
  A1(1,0) = 3.0; A1(1,1) = 2.0;

  B1(0,0) = 6.0;
  B1(1,0) = 8.0;

  Lucee::accumulate(C1, A1, B1);
// check if accumulate worked
  LC_ASSERT("Testing if accumulate worked", C1(0,0) == 20.0);
  LC_ASSERT("Testing if accumulate worked", C1(1,0) == 34.0);
}

void
test_8()
{
  Lucee::Matrix<double> A(2,3);
  Lucee::Vector<double> x(3), y(2);

  A = 2.0;
  x = 1.0;
  y = 0.0;

  accumulate(y, A, x);
  LC_ASSERT("Testing matrix-vector product", y[0] == 6.0);
  LC_ASSERT("Testing matrix-vector product", y[1] == 6.0);

  Lucee::Matrix<double> A1(2,3);
  Lucee::Vector<double> x1(3), y1(2);

  A1 = 2.0;
  x1[0] = 1.0; x1[1] = 2.0; x1[2] = 3.0;
  accumulate(y1, A1, x1);
  LC_ASSERT("Testing matrix-vector product", y1[0] == 12.0);
  LC_ASSERT("Testing matrix-vector product", y1[1] == 12.0);
}

void
test_9()
{
  Lucee::Vector<double> x(3), y(2);
  Lucee::Matrix<double> A(3, 2);

  A = 0.0;
  x = 1.0;
  y = 2.0;
  accumulate(A, 1.0, x, y);
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing if vector-vector outer product worked", A(i,j)==2.0);
}

void
test_10()
{
  Lucee::Matrix<double> A(3, 5);
  A = 10.5;
  for (int i=0; i<5; ++i)
    A(2,i) = (i+0.5)*2.5;

  LC_RAISES("Testing if exception if thrown", A.getRow(3), Lucee::Except);
  LC_RAISES("Testing if exception if thrown", A.getCol(5), Lucee::Except);

  Lucee::Vector<double> x = A.getRow(2);
  LC_ASSERT("Testing if row of matrix got correctly", x.getLength() == 5);

  for (int i=0; i<5; ++i)
    LC_ASSERT("Testing if row of matrix has proper data", x(i) == (i+0.5)*2.5);

  x = 12.5;
  for (int i=0; i<5; ++i)
    LC_ASSERT("Testing if matrix changed", A(2,i) == 12.5);

  Lucee::Vector<double> y = A.getCol(4);
  LC_ASSERT("Testing if row of matrix got correctly", y.getLength() == 3);

  for (int i=0; i<3; ++i)
    A(i,4) = (i+0.5)*2.5;

  for (int i=0; i<3; ++i)
    LC_ASSERT("Testing if row of matrix has proper data", y(i) == (i+0.5)*2.5);

  y = 12.5;
  for (int i=0; i<3; ++i)
    LC_ASSERT("Testing if matrix changed", A(i,4) == 12.5);
}

void
test_11()
{
  Lucee::Matrix<double> A(3,3);
  A = 0.0; // everything is zero, but ...
  A(0,0) = 3;
  A(1,1) = 2;
  A(2,2) = 12.5;

  Lucee::Matrix<double> B(3,1);
  B(0,0) = 6.0;
  B(1,0) = 4.0;
  B(2,0) = 25.0;

// solve linear system
  Lucee::solve(A, B);
// test solution
  LC_ASSERT("Testing solution to linear system", B(0,0) == 2);
  LC_ASSERT("Testing solution to linear system", B(1,0) == 2);
  LC_ASSERT("Testing solution to linear system", B(2,0) == 2);
}

void
test_12()
{
  Lucee::Matrix<double> A(10, 12);
  Lucee::Vector<double> rowFact(10), colFact(12);

  for (int i=rowFact.getLower(0); i<rowFact.getUpper(0); ++i)
    rowFact[i] = i;

  for (int i=colFact.getLower(0); i<colFact.getUpper(0); ++i)
    colFact[i] = i;

  A = 1.0;
  A.scaleRows(rowFact);
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing if row-scaling worked", A(i,j) == rowFact[i]);

  A = 1.0;
  A.scaleCols(colFact);
  for (int i=A.getLower(0); i<A.getUpper(0); ++i)
    for (int j=A.getLower(1); j<A.getUpper(1); ++j)
      LC_ASSERT("Testing if col-scaling worked", A(i,j) == colFact[j]);
}

void
test_13()
{
  Lucee::Matrix<double> A(10, 10);
  A = 2.0;

// get view
  Lucee::Matrix<double> view = A.getView(0, 5, 0, 5);

  LC_ASSERT("Testing view bounds", view.getLower(0) == 0);
  LC_ASSERT("Testing view bounds", view.getUpper(0) == 5);

  LC_ASSERT("Testing view bounds", view.getLower(1) == 0);
  LC_ASSERT("Testing view bounds", view.getUpper(1) == 5);

  for (int i=view.getLower(0); i<view.getUpper(0); ++i)
    for (int j=view.getLower(1); j<view.getUpper(1); ++j)
      view(i,j) = (10*i+j)*0.5;

  for (int i=view.getLower(0); i<view.getUpper(0); ++i)
    for (int j=view.getLower(1); j<view.getUpper(1); ++j)
      LC_ASSERT("Testing view values", A(i,j) == (10*i+j)*0.5);

// get view
  Lucee::Matrix<double> view2 = A.getView(5, 10, 0, 5);

  LC_ASSERT("Testing view bounds", view.getLower(0) == 0);
  LC_ASSERT("Testing view bounds", view.getUpper(0) == 5);

  LC_ASSERT("Testing view bounds", view.getLower(1) == 0);
  LC_ASSERT("Testing view bounds", view.getUpper(1) == 5);

  for (int i=view.getLower(0); i<view.getUpper(0); ++i)
    for (int j=view.getLower(1); j<view.getUpper(1); ++j)
      view(i,j) = (10*i+j)*0.5;

  for (int i=view.getLower(0); i<view.getUpper(0); ++i)
    for (int j=view.getLower(1); j<view.getUpper(1); ++j)
      LC_ASSERT("Testing view values", A(i,j) == (10*i+j)*0.5);
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
  test_6();
  test_7();
  test_8();
  test_9();
  test_10();
  test_11();
  test_12();
  test_13();
  LC_END_TESTS;
}
