/**
 * @file	lcfixedvector.cxx
 *
 * @brief	Unit tests for Lucee::FixedVector class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcFixedVector.h>
#include <LcTest.h>

// std includes
#include <cmath>

double
computeNorm(unsigned n, double *vec)
{
  double sum = 0.0;
  for (unsigned i=0; i<n; ++i)
    sum += vec[i]*vec[i];
  return sqrt(sum);
}

template <unsigned NELEM>
double
computeNorm2(double vec[NELEM])
{
  double sum = 0.0;
  for (unsigned i=0; i<NELEM; ++i)
    sum += vec[i]*vec[i];
  return sqrt(sum);
}

void
test_1()
{
  Lucee::FixedVector<3, double> xv(1.0);

  LC_ASSERT("Testing size of vector", xv.numElem() == 3);

  for (unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing value of fixed vector", xv[i]==1.0);

  for (unsigned i=0; i<3; ++i)
    xv[i] = i;

  for (unsigned i=0; i<3; ++i)
    LC_ASSERT("Testing value of fixed vector", xv[i]==i);

  double sum = 0.0;
  for (unsigned i=0; i<3; ++i)
    sum += xv[i]*xv[i];
  double myNorm = sqrt(sum);

// this tests if we can call method which takes a double*
  double norm = computeNorm(xv.numElem(), &xv[0]);
// this tests if we can call method which takes a double[NELEM]
  double norm2 = computeNorm2<3>(&xv[0]);

  LC_ASSERT("Testing norm", norm==myNorm);
  LC_ASSERT("Testing norm", norm2==myNorm);
}

void
test_2()
{
  Lucee::FixedVector<3, double> xv(1.0, 2.0, 3.0);

  LC_ASSERT("Testing first entry in array", xv[0]==1.0);
  LC_ASSERT("Testing second entry in array", xv[1]==2.0);
  LC_ASSERT("Testing third entry in array", xv[2]==3.0);

  double norm = sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0);
  LC_ASSERT("Testing norm of vector", epsCmp(xv.getNorm(), norm));
}

int
main(void)
{
  LC_BEGIN_TESTS("lcfixedvector");
  test_1();
  test_2();
  LC_END_TESTS;
}
