/**
 * @file	lcarray.cxx
 *
 * @brief	Unit tests for Lucee::Array class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcArray.h>
#include <LcTest.h>

void
test_1()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);

  int start[2] = {2, 0};
  Lucee::Array<2, double> arr2(shape, start);

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking start indices", arr1.getLower(i) == 0);
    LC_ASSERT("Checking start indices", arr2.getLower(i) == start[i]);

    LC_ASSERT("Checking end indices", arr1.getUpper(i) == (int) shape[i]);
    int end = start[i]+shape[i];
    LC_ASSERT("Checking end indices", arr2.getUpper(i) == end);

    LC_ASSERT("Checking shape", arr1.getShape(i) == shape[i]);
  }

  unsigned myShape[2];
  arr2.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  arr1.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  LC_ASSERT("Testing size of array", 50 == arr1.getSize());
  LC_ASSERT("Testing size of array", 50 == arr2.getSize());

  LC_ASSERT("Testing if array is contigous", true == arr1.isContiguous());
  LC_ASSERT("Testing if array is contigous", true == arr2.isContiguous());
}

void
test_1d()
{
  unsigned shape[1] = {10};
  Lucee::Array<1, double> arr1(shape);

  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    arr1(i) = 10*i;

  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", arr1(i) == 10*i);

  int idx[1];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    {
      idx[0] = i;
      LC_ASSERT("Testing 2D indexer", arr1(idx) == 10*i);
    }
}

void
test_2d()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);

  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      arr1(i,j) = 10*i+j;

  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      LC_ASSERT("Testing 2D indexer", arr1(i,j) == 10*i+j);

  int idx[2];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", arr1(idx) == 10*i+j);
    }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcarray");
  test_1();
  test_1d();
  test_2d();
  LC_END_TESTS;
}
