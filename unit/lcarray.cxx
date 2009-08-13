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

  double count = 0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    arr1(i) = count++;

  count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", arr1(i) == count++);

  count = 0.0;
  int idx[1];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", arr1(idx) == count++);
  }
}

void
test_2d()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);

  double count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      arr1(i,j) = count++;

  count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
    {
      LC_ASSERT("Testing 2D indexer", arr1(i,j) == count++);
    }

  count = 0.0;
  int idx[2];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", arr1(idx) == count++);
    }
}

void
test_3d()
{
  unsigned shape[3] = {5, 10, 12};
  Lucee::Array<3, double> arr1(shape);

  double count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
        arr1(i,j,k) = count++;

  count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
      {
        LC_ASSERT("Testing 3D indexer", arr1(i,j,k) == count++);
      }

  count = 0.0;
  int idx[3];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indexer", arr1(idx) == count++);
      }
}

void
test_4d()
{
  unsigned shape[4] = {5, 10, 12, 6};
  Lucee::Array<4, double> arr1(shape);

  double count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
        for (int l=arr1.getLower(3); l<arr1.getUpper(3); ++l)
          arr1(i,j,k,l) = count++;

  count = 0.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
        for (int l=arr1.getLower(3); l<arr1.getUpper(3); ++l)
        {
          LC_ASSERT("Testing 4D indexer", arr1(i,j,k,l) == count++);
        }

  count = 0.0;
  int idx[4];
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      for (int k=arr1.getLower(2); k<arr1.getUpper(2); ++k)
        for (int l=arr1.getLower(3); l<arr1.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", arr1(idx) == count++);
        }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcarray");
  test_1();
  test_1d();
  test_2d();
  test_3d();
  test_4d();
  LC_END_TESTS;
}
