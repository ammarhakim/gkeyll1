/**
 * @file	lcarray.cxx
 *
 * @brief	Unit tests for Lucee::Array class
 */

// lucee includes
#include <LcArray.h>
#include <LcRowMajorSequencer.h>
#include <LcTest.h>

// std includes
#include <complex>

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

  int myStart[2];
  arr2.fillWithStart(myStart);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking start", myStart[i] == start[i]);

  arr1.fillWithStart(myStart);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myStart[i] == 0);

  LC_ASSERT("Testing size of array", 50 == arr1.getSize());
  LC_ASSERT("Testing size of array", 50 == arr2.getSize());

  LC_ASSERT("Testing if array is contigous", true == arr1.isContiguous());
  LC_ASSERT("Testing if array is contigous", true == arr2.isContiguous());

  arr1 = 10.0;
  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arr1(i,j)==10.0);

  arr2 = 10.0;
  for (int i=arr2.getLower(0); i<arr2.getUpper(0); ++i)
    for (int j=arr2.getLower(1); j<arr2.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arr2(i,j)==10.0);

  Lucee::Region<2, int> rgn = arr2.getRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing region returned by array",
      rgn.getLower(i) == start[i]);
    LC_ASSERT("Testing region returned by array",
      rgn.getUpper(i) == (start[i]+shape[i]));
  }
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

void
test_5()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);
  int start[2] = {2, 0};
  Lucee::Array<2, double> arr2(shape, start);

// call copy ctor
  Lucee::Array<2, double> arrCopy(arr1);
  Lucee::Array<2, double> arr2Copy(arr2);

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking start indices", arrCopy.getLower(i) == 0);
    LC_ASSERT("Checking start indices", arr2Copy.getLower(i) == start[i]);

    LC_ASSERT("Checking end indices", arrCopy.getUpper(i) == (int) shape[i]);
    int end = start[i]+shape[i];
    LC_ASSERT("Checking end indices", arr2Copy.getUpper(i) == end);

    LC_ASSERT("Checking shape", arrCopy.getShape(i) == shape[i]);
  }

  unsigned myShape[2];
  arr2Copy.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  arrCopy.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  LC_ASSERT("Testing size of array", 50 == arrCopy.getSize());
  LC_ASSERT("Testing size of array", 50 == arr2Copy.getSize());

  LC_ASSERT("Testing if array is contigous", true == arrCopy.isContiguous());
  LC_ASSERT("Testing if array is contigous", true == arr2Copy.isContiguous());

  arrCopy = 10.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arrCopy(i,j)==10.0);

  arr2Copy = 10.0;
  for (int i=arr2Copy.getLower(0); i<arr2Copy.getUpper(0); ++i)
    for (int j=arr2Copy.getLower(1); j<arr2Copy.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arr2Copy(i,j)==10.0);

  double count = 0.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
      arrCopy(i,j) = count++;

  count = 0.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
    {
      LC_ASSERT("Testing 2D indexer", arrCopy(i,j) == count++);
    }

  count = 0.0;
  int idx[2];
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", arrCopy(idx) == count++);
    }
}

void
test_6()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr(shape);
  arr = 10.0;

  Lucee::Array<2, double> arrCopy(arr);

  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
      LC_ASSERT("Testing if original affects copy", arrCopy(i,j)==10.0);

  arrCopy = 20.0;

  for (int i=arr.getLower(0); i<arr.getUpper(0); ++i)
    for (int j=arr.getLower(1); j<arr.getUpper(1); ++j)
      LC_ASSERT("Testing if original affects copy", arr(i,j)==20.0);

  arrCopy(2,2) = 12.5;
  LC_ASSERT("Testing if copy is shallow copy", arr(2,2) == 12.5);
}

void
test_7()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);
  int start[2] = {2, 0};
  Lucee::Array<2, double> arr2(shape, start);

// create small 2D arrays
  int zeros[2] = {0, 0};
  unsigned ones[2] = {1, 1};
  Lucee::Array<2, double> arrCopy(ones, zeros);
  Lucee::Array<2, double> arr2Copy(ones, zeros);

// test assignment operator
  arrCopy = arr1;
  arr2Copy = arr2;

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking start indices", arrCopy.getLower(i) == 0);
    LC_ASSERT("Checking start indices", arr2Copy.getLower(i) == start[i]);

    LC_ASSERT("Checking end indices", arrCopy.getUpper(i) == (int) shape[i]);
    int end = start[i]+shape[i];
    LC_ASSERT("Checking end indices", arr2Copy.getUpper(i) == end);

    LC_ASSERT("Checking shape", arrCopy.getShape(i) == shape[i]);
  }

  unsigned myShape[2];
  arr2Copy.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  arrCopy.fillWithShape(myShape);
  for (unsigned i=0; i<2; ++i)
    LC_ASSERT("Checking shape", myShape[i] == shape[i]);

  LC_ASSERT("Testing size of array", 50 == arrCopy.getSize());
  LC_ASSERT("Testing size of array", 50 == arr2Copy.getSize());

  LC_ASSERT("Testing if array is contigous", true == arrCopy.isContiguous());
  LC_ASSERT("Testing if array is contigous", true == arr2Copy.isContiguous());

  arrCopy = 10.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arrCopy(i,j)==10.0);

  arr2Copy = 10.0;
  for (int i=arr2Copy.getLower(0); i<arr2Copy.getUpper(0); ++i)
    for (int j=arr2Copy.getLower(1); j<arr2Copy.getUpper(1); ++j)
      LC_ASSERT("Testing is setting all values worked", arr2Copy(i,j)==10.0);

  double count = 0.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
      arrCopy(i,j) = count++;

  count = 0.0;
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
    {
      LC_ASSERT("Testing 2D indexer", arrCopy(i,j) == count++);
    }

  count = 0.0;
  int idx[2];
  for (int i=arrCopy.getLower(0); i<arrCopy.getUpper(0); ++i)
    for (int j=arrCopy.getLower(1); j<arrCopy.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", arrCopy(idx) == count++);
    }
}

void
test_8()
{
  int lower[4] = {0, 0, 0};
  int upper[4] = {5, 10, 12, 6};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::Array<4, double> arr1(rgn);

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

void
test_9()
{
  int lower[2] = {0, 0};
  int upper[2] = {15, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Array<2, double> arr(rgn, 12.5);

  int loSl[2] = {2, 3};
  int upSl[2] = {12, 8};
  Lucee::Region<2, int> rgnSl(loSl, upSl);
// slice it
  Lucee::Array<2, double> arrSl = arr.getSlice(rgnSl);
// check slice
  LC_ASSERT("Testing slice lower bounds", arrSl.getLower(0) == 2);
  LC_ASSERT("Testing slice lower bounds", arrSl.getLower(1) == 3);

  LC_ASSERT("Testing slice upper bounds", arrSl.getUpper(0) == 12);
  LC_ASSERT("Testing slice upper bounds", arrSl.getUpper(1) == 8);

// set complete sliced array
  for (int i=arrSl.getLower(0); i<arrSl.getUpper(0); ++i)
    for (int j=arrSl.getLower(1); j<arrSl.getUpper(1); ++j)
      arrSl(i,j) = 18.5;

  for (int i=arrSl.getLower(0); i<arrSl.getUpper(0); ++i)
    for (int j=arrSl.getLower(1); j<arrSl.getUpper(1); ++j)
      LC_ASSERT("Testing if setting slice worked", arr(i,j) == 18.5);
}

void
test_10()
{
  int lower[2] = {0, 0};
  int upper[2] = {15, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Array<2, double> arr(rgn, 12.5);

  int loSl[2] = {2, 3};
  int upSl[2] = {12, 8};
  Lucee::Region<2, int> rgnSl(loSl, upSl);
// slice it
  Lucee::Array<2, double> arrSl = arr.getSlice(rgnSl);
// check slice
  LC_ASSERT("Testing slice lower bounds", arrSl.getLower(0) == 2);
  LC_ASSERT("Testing slice lower bounds", arrSl.getLower(1) == 3);

  LC_ASSERT("Testing slice upper bounds", arrSl.getUpper(0) == 12);
  LC_ASSERT("Testing slice upper bounds", arrSl.getUpper(1) == 8);

// set complete sliced array
  arrSl = 18.5;
  for (int i=arrSl.getLower(0); i<arrSl.getUpper(0); ++i)
    for (int j=arrSl.getLower(1); j<arrSl.getUpper(1); ++j)
      arrSl(i,j) = 18.5;
}

void
test_11()
{
  int lower[2] = {0, 0};
  int upper[2] = {15, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Array<2, double> arr(rgn, 12.5);

  int loSl[2] = {2, 3};
  int upSl[2] = {17, 8};
  Lucee::Region<2, int> rgnSl(loSl, upSl);
// try to slice it
  LC_RAISES("Testing if slicing invalid region works", arr.getSlice(rgnSl), Lucee::Except);
}

void
test_12()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::Array<4, double> arr(rgn, 12.5);

// deflate the 4D array into a 2D array
  unsigned defDims[2] = {1, 2};
  int defDimsIdx[2] = {10, 10};

  Lucee::Array<2, double> defArr = arr.deflate<2>(defDims, defDimsIdx);

  LC_ASSERT("Testing deflated array lower bound", defArr.getLower(0) == 3);
  LC_ASSERT("Testing deflated array lower bound", defArr.getLower(1) == 10);

  LC_ASSERT("Testing deflated array upper bound", defArr.getUpper(0) == 15);
  LC_ASSERT("Testing deflated array upper bound", defArr.getUpper(1) == 30);

// set the 4D array
  Lucee::RowMajorSequencer<4> seq(rgn);
  double count = 1.5;
  while (seq.step())
  {
    arr(seq.getIndex()) = count;
    count += 1.5;
  }

// now check deflated array
  for (int i=defArr.getLower(0); i<defArr.getUpper(0); ++i)
    for (int j=defArr.getLower(1); j<defArr.getUpper(1); ++j)
      LC_ASSERT("Testing deflated array", defArr(i,j) == arr(i,10,10,j));
}

void
test_13()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::Array<4, double> arr(rgn, 12.5);
// make duplicate
  Lucee::Array<4, double> arrDup = arr.duplicate();
  
// test it
  for (unsigned i=0; i<4; ++i)
  {
    LC_ASSERT("Testing duplicate array", arrDup.getLower(i) == arr.getLower(i));
    LC_ASSERT("Testing duplicate array", arrDup.getUpper(i) == arr.getUpper(i));
    LC_ASSERT("Testing duplicate array", arrDup.getShape(i) == arr.getShape(i));
  }

  Lucee::ColMajorSequencer<4> seq(rgn);
  while (seq.step())
    LC_ASSERT("Testing duplicate array", arrDup(seq.getIndex()) == arr(seq.getIndex()));
}

void
test_14()
{
  int lower[1] = {3};
  int upper[1] = {3};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::Array<1, double> zeroVol(rgn);

  LC_ASSERT("Testing zero-volume region", rgn.getVolume() == 0);
}

void
test_15()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);

  int start[2] = {2, 0};
  Lucee::Array<2, double> arr2(shape, start);

  for (int i=arr1.getLower(0); i<arr1.getUpper(0); ++i)
    for (int j=arr1.getLower(1); j<arr1.getUpper(1); ++j)
      arr1(i,j) = (5*i+j)*0.5;

// copy to arr2
  arr2.copy(arr1);
  for (int i=arr1.getLower(0), i1=arr2.getLower(0); i<arr1.getUpper(0); ++i, ++i1)
    for (int j=arr1.getLower(1), j1=arr2.getLower(1); j<arr1.getUpper(1); ++j, ++j1)
      LC_ASSERT("Testing if array copy worked", arr1(i,j) == arr2(i1,j1));

  arr2 = 0.0;
  arr2.accumulate(0.5, arr1);
  for (int i=arr1.getLower(0), i1=arr2.getLower(0); i<arr1.getUpper(0); ++i, ++i1)
    for (int j=arr1.getLower(1), j1=arr2.getLower(1); j<arr1.getUpper(1); ++j, ++j1)
      LC_ASSERT("Testing if array accumulate worked", 0.5*arr1(i,j) == arr2(i1,j1));

  arr2.accumulate(1.0, arr2);
  for (int i=arr1.getLower(0), i1=arr2.getLower(0); i<arr1.getUpper(0); ++i, ++i1)
    for (int j=arr1.getLower(1), j1=arr2.getLower(1); j<arr1.getUpper(1); ++j, ++j1)
      LC_ASSERT("Testing if array accumulate worked", arr1(i,j) == arr2(i1,j1));
}

void
test_16()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr(shape);

  arr = 22.5;
  arr += 5.25;
  for (int i=arr.getLower(0); i<arr.getUpper(0); ++i)
    for (int j=arr.getLower(1); j<arr.getUpper(1); ++j)
      LC_ASSERT("Testing if += worked correctly", arr(i,j) == 22.5+5.25);

  arr = 22.5;
  arr *= 5.25;
  for (int i=arr.getLower(0); i<arr.getUpper(0); ++i)
    for (int j=arr.getLower(1); j<arr.getUpper(1); ++j)
      LC_ASSERT("Testing if *= worked correctly", arr(i,j) == 22.5*5.25);

  arr = 22.5;
  arr /= 5.25;
  for (int i=arr.getLower(0); i<arr.getUpper(0); ++i)
    for (int j=arr.getLower(1); j<arr.getUpper(1); ++j)
      LC_ASSERT("Testing if += worked correctly", arr(i,j) == 22.5/5.25);
}

void
test_17()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, std::complex<double> > arr(shape);
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcarray");
  test_1();
  test_1d();
  test_2d();
  test_3d();
  test_4d();
  test_5();
  test_6();
  test_7();
  test_8();
  test_9();
  test_10();
  test_11();
  test_12();
  test_13();
  test_14();
  test_15();
  test_16();
  test_17();
  LC_END_TESTS;
}
