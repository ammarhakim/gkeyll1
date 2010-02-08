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
test_8()
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
  LC_END_TESTS;
}
