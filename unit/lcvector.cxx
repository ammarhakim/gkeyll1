/**
 * @file	lcvector.cxx
 *
 * @brief	Unit tests for Lucee::Vector class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcArray.h>
#include <LcVector.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  Lucee::Vector<double> vec(10);

  LC_ASSERT("Testing size of vector", vec.getLength() == 10);
  vec = 10.5;

  for (unsigned i=0; i<vec.getLength(); ++i)
    LC_ASSERT("Testing value of fixed vector", vec[i] == 10.5);

  for (unsigned i=0; i<vec.getLength(); ++i)
    vec[i] = i;

  for (unsigned i=0; i<vec.getLength(); ++i)
    LC_ASSERT("Testing value of fixed vector", vec[i] == i);
}

void
test_2()
{
  int lower[1] = {10};
  int upper[1] = {20};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::Array<1, double> arr1d(rgn);

  for (int i=arr1d.getLower(0); i<arr1d.getUpper(0); ++i)
    arr1d(i) = (0.5+i)*2.5;

// create vector from array
  Lucee::Vector<double> vec(arr1d);
  for (int i=arr1d.getLower(0); i<arr1d.getUpper(0); ++i)
    LC_ASSERT("Testing vector created from 1D array", vec(i) == (0.5+i)*2.5);
}

void
test_3()
{
  Lucee::Vector<double> vec(20);
  vec = 25.0;
  vec.scale(0.5);
  for (int i=vec.getLower(0); i<vec.getUpper(0); ++i)
    LC_ASSERT("Testing scale() function", vec[i] == 12.5);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcvector");
  test_1();
  test_2();
  test_3();
  LC_END_TESTS;
}
