/**
 * @file	lccolmajorindexer.cxx
 *
 * @brief	Unit tests for Lucee::Array class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcColMajorIndexer.h>
#include <LcTest.h>

void
test_1()
{
  int start[1] = {1};
  unsigned shape[1] = {10};
  Lucee::ColMajorIndexer<1> col(shape, start);

  int count=0;
  for (int i=col.getLower(0); i<col.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", col.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=col.getLower(0); i<col.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", col.getGenIndex(idx) == count++);
  }
}

void
test_2()
{
  int start[2] = {1, 2};
  unsigned shape[2] = {10, 15};
  Lucee::ColMajorIndexer<2> col(shape, start);

  int count=0;
  for (int j=col.getLower(1); j<col.getUpper(1); ++j)
    for (int i=col.getLower(0); i<col.getUpper(0); ++i)
      LC_ASSERT("Testing 2D indexer", col.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int j=col.getLower(1); j<col.getUpper(1); ++j)
    for (int i=col.getLower(0); i<col.getUpper(0); ++i)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", col.getGenIndex(idx) == count++);
    }
}

void
test_3()
{
  int start[3] = {1, 2, -2};
  unsigned shape[3] = {10, 15, 20};
  Lucee::ColMajorIndexer<3> col(shape, start);

  int count=0;
  for (int k=col.getLower(2); k<col.getUpper(2); ++k)
    for (int j=col.getLower(1); j<col.getUpper(1); ++j)
      for (int i=col.getLower(0); i<col.getUpper(0); ++i)
        LC_ASSERT("Testing 3D indexer", col.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];
  for (int k=col.getLower(2); k<col.getUpper(2); ++k)
    for (int j=col.getLower(1); j<col.getUpper(1); ++j)
      for (int i=col.getLower(0); i<col.getUpper(0); ++i)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indexer", col.getGenIndex(idx) == count++);
      }
}

void
test_4()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::ColMajorIndexer<4> col(shape, start);

  int count=0;
  for (int l=col.getLower(3); l<col.getUpper(3); ++l)
    for (int k=col.getLower(2); k<col.getUpper(2); ++k)
      for (int j=col.getLower(1); j<col.getUpper(1); ++j)
        for (int i=col.getLower(0); i<col.getUpper(0); ++i)
          LC_ASSERT("Testing 4D indexer", col.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int l=col.getLower(3); l<col.getUpper(3); ++l)
    for (int k=col.getLower(2); k<col.getUpper(2); ++k)
      for (int j=col.getLower(1); j<col.getUpper(1); ++j)
        for (int i=col.getLower(0); i<col.getUpper(0); ++i)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", col.getGenIndex(idx) == count++);
        }
}

void
test_5()
{
  int start[2] = {1, 2};
  unsigned shape[2] = {10, 15};
  Lucee::ColMajorIndexer<2> col(shape, start);
  Lucee::ColMajorIndexer<2> colCopy(col); // copy it
 
// now test the copied indexer (testing copy ctor)
  int count=0;
  for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
    for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
      LC_ASSERT("Testing 2D indexer", colCopy.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
    for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", colCopy.getGenIndex(idx) == count++);
    }
}

void
test_6()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::ColMajorIndexer<4> col(shape, start);
  Lucee::ColMajorIndexer<4> colCopy(col);

  int count=0;
  for (int l=colCopy.getLower(3); l<colCopy.getUpper(3); ++l)
    for (int k=colCopy.getLower(2); k<colCopy.getUpper(2); ++k)
      for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
        for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
          LC_ASSERT("Testing 4D indexer", colCopy.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int l=colCopy.getLower(3); l<colCopy.getUpper(3); ++l)
    for (int k=colCopy.getLower(2); k<colCopy.getUpper(2); ++k)
      for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
        for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", colCopy.getGenIndex(idx) == count++);
        }
}

void
test_7()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::ColMajorIndexer<4> col(shape, start);

  int zeros[4] = {0, 0, 0, 0};
  unsigned ones[4] = {1, 1, 1, 1};
  Lucee::ColMajorIndexer<4> colCopy(ones, zeros);
  colCopy = col;

  int count=0;
  for (int l=colCopy.getLower(3); l<colCopy.getUpper(3); ++l)
    for (int k=colCopy.getLower(2); k<colCopy.getUpper(2); ++k)
      for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
        for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
          LC_ASSERT("Testing 4D indexer", colCopy.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int l=colCopy.getLower(3); l<colCopy.getUpper(3); ++l)
    for (int k=colCopy.getLower(2); k<colCopy.getUpper(2); ++k)
      for (int j=colCopy.getLower(1); j<colCopy.getUpper(1); ++j)
        for (int i=colCopy.getLower(0); i<colCopy.getUpper(0); ++i)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", colCopy.getGenIndex(idx) == count++);
        }
}

void
test_1_1()
{
  int lower[1] = {1};
  int upper[1] = {10+1};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<1> col(rgn);

  int count=0;
  for (int i=col.getLower(0); i<col.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", col.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=col.getLower(0); i<col.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", col.getGenIndex(idx) == count++);
  }
}

void
test_2_2()
{
  int lower[2] = {1, 2};
  int upper[2] = {10+1, 15+2};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<2> col(rgn);

  int count=0;
  for (int j=col.getLower(1); j<col.getUpper(1); ++j)
    for (int i=col.getLower(0); i<col.getUpper(0); ++i)
      LC_ASSERT("Testing 2D indexer", col.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int j=col.getLower(1); j<col.getUpper(1); ++j)
    for (int i=col.getLower(0); i<col.getUpper(0); ++i)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", col.getGenIndex(idx) == count++);
    }
}

void
test_3_3()
{
  int lower[3] = {1, 2, -2};
  int upper[3] = {10+1, 15+2, 20-2};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<3> col(rgn);

  int count=0;
  for (int k=col.getLower(2); k<col.getUpper(2); ++k)
    for (int j=col.getLower(1); j<col.getUpper(1); ++j)
      for (int i=col.getLower(0); i<col.getUpper(0); ++i)
        LC_ASSERT("Testing 3D indexer", col.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];
  for (int k=col.getLower(2); k<col.getUpper(2); ++k)
    for (int j=col.getLower(1); j<col.getUpper(1); ++j)
      for (int i=col.getLower(0); i<col.getUpper(0); ++i)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indexer", col.getGenIndex(idx) == count++);
      }
}

void
test_4_4()
{
  int lower[4] = {1, 2, -2, 3};
  int upper[4] = {10+1, 15+2, 20-2, 12+3};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<4> col(rgn);

  int count=0;
  for (int l=col.getLower(3); l<col.getUpper(3); ++l)
    for (int k=col.getLower(2); k<col.getUpper(2); ++k)
      for (int j=col.getLower(1); j<col.getUpper(1); ++j)
        for (int i=col.getLower(0); i<col.getUpper(0); ++i)
          LC_ASSERT("Testing 4D indexer", col.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int l=col.getLower(3); l<col.getUpper(3); ++l)
    for (int k=col.getLower(2); k<col.getUpper(2); ++k)
      for (int j=col.getLower(1); j<col.getUpper(1); ++j)
        for (int i=col.getLower(0); i<col.getUpper(0); ++i)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", col.getGenIndex(idx) == count++);
        }
}

template <unsigned NDIM>
void
createRowMajorIndexer(const unsigned shape[NDIM], const int start[NDIM], int ai[NDIM+1])
{
  ai[NDIM] = 1;
  for (unsigned i=NDIM-1; i>=1; --i)
    ai[i] = ai[i+1]*shape[i];
  
  int sum = 0;
  for (unsigned i=1; i<NDIM+1; ++i)
    sum += ai[i]*start[i-1];
  ai[0] = -sum;
}

void
test_1_r()
{
  int start[1] = {1};
  unsigned shape[1] = {10};
  int ai[2];
  createRowMajorIndexer<1>(shape, start, ai);
// fake row major indexing
  Lucee::LinIndexer<1> lin(shape, start, ai);

  int count=0;
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", lin.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", lin.getGenIndex(idx) == count++);
  }
}

void
test_2_r()
{
  int start[2] = {1, 2};
  unsigned shape[2] = {10, 15};
  int ai[3];
  createRowMajorIndexer<2>(shape, start, ai);
  
  Lucee::LinIndexer<2> lin(shape, start, ai);
  int count=0;
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
      LC_ASSERT("Testing 2D row-lin indexer", lin.getIndex(i,j) == count++);
  
  count=0;
 int idx[2];
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D row-lin indexer", lin.getGenIndex(idx) == count++);
    }
}

void
test_3_r()
{
  int start[3] = {1, 2, -2};
  unsigned shape[3] = {10, 15, 20};
  int ai[4];
  createRowMajorIndexer<3>(shape, start, ai);

  Lucee::LinIndexer<3> lin(shape, start, ai);

  int count=0;
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
      for (int k=lin.getLower(2); k<lin.getUpper(2); ++k)
        LC_ASSERT("Testing 3D row-lin indexer", lin.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];

  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
      for (int k=lin.getLower(2); k<lin.getUpper(2); ++k)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D row-lin indexer", lin.getGenIndex(idx) == count++);
      }
}

void
test_4_r()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  int ai[5];
  createRowMajorIndexer<4>(shape, start, ai);

  Lucee::LinIndexer<4> lin(shape, start, ai);

  int count=0;
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
      for (int k=lin.getLower(2); k<lin.getUpper(2); ++k)
        for (int l=lin.getLower(3); l<lin.getUpper(3); ++l)
          LC_ASSERT("Testing 4D row-lin indexer", lin.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int i=lin.getLower(0); i<lin.getUpper(0); ++i)
    for (int j=lin.getLower(1); j<lin.getUpper(1); ++j)
      for (int k=lin.getLower(2); k<lin.getUpper(2); ++k)
        for (int l=lin.getLower(3); l<lin.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D row-lin indexer", lin.getGenIndex(idx) == count++);
        }
}

void
test_8()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<4> col(rgn);

// deflate the 4D indexer into a 2D indexer
  unsigned defDims[2] = {1, 2};
  int defDimsIdx[2] = {10, 10};

  Lucee::ColMajorIndexer<2> defCol
    = col.deflate<2>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defCol.getLower(0) == 3);
  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defCol.getLower(1) == 10);

  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defCol.getUpper(0) == 15);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defCol.getUpper(1) == 30);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      LC_ASSERT("Testing deflated indexer", defCol.getIndex(i,j) == col.getIndex(i,10,10,j));
}

void
test_9()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<4> col(rgn);

// deflate the 4D indexer into a 2D indexer
  unsigned defDims[2] = {0, 3};
  int defDimsIdx[2] = {4, 15};

  Lucee::ColMajorIndexer<2> defCol
    = col.deflate<2>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defCol.getLower(0) == 6);
  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defCol.getLower(1) == 8);

  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defCol.getUpper(0) == 20);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defCol.getUpper(1) == 25);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      LC_ASSERT("Testing deflated indexer", defCol.getIndex(i,j) == col.getIndex(4,i,j,15));
}

void
test_10()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::ColMajorIndexer<4> col(rgn);

// deflate the 4D indexer into a 1D indexer
  unsigned defDims[3] = {0, 1, 2};
  int defDimsIdx[3] = {4, 15, 20};

  Lucee::ColMajorIndexer<1> defCol
    = col.deflate<1>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defCol.getLower(0) == 10);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defCol.getUpper(0) == 30);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    LC_ASSERT("Testing deflated indexer", defCol.getIndex(i) == col.getIndex(4,15,20,i));
}

int
main(void) 
{
  LC_BEGIN_TESTS("lccolmajorindexer");
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

  test_1_1();
  test_2_2();
  test_3_3();
  test_4_4();

  test_1_r();
  test_2_r();
  test_3_r();
  test_4_r();
  LC_END_TESTS;
}
