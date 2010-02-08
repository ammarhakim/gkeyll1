/**
 * @file	lcrowmajorindexer.cxx
 *
 * @brief	Unit tests for Lucee::Array class
 *
 * @version	$Id: lcrowmajorindexer.cxx 142 2009-08-23 16:32:29Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcRowMajorIndexer.h>
#include <LcTest.h>

void
test_1()
{
  int start[1] = {1};
  unsigned shape[1] = {10};
  Lucee::RowMajorIndexer<1> row(shape, start);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", row.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", row.getGenIndex(idx) == count++);
  }
}

void
test_2()
{
  int start[2] = {1, 2};
  unsigned shape[2] = {10, 15};
  Lucee::RowMajorIndexer<2> row(shape, start);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
  for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      LC_ASSERT("Testing 2D indexer", row.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", row.getGenIndex(idx) == count++);
    }
}

void
test_3()
{
  int start[3] = {1, 2, -2};
  unsigned shape[3] = {10, 15, 20};
  Lucee::RowMajorIndexer<3> row(shape, start);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        LC_ASSERT("Testing 3D indexer", row.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];

  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indexer", row.getGenIndex(idx) == count++);
      }
}

void
test_4()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::RowMajorIndexer<4> row(shape, start);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        for (int l=row.getLower(3); l<row.getUpper(3); ++l)
          LC_ASSERT("Testing 4D indexer", row.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        for (int l=row.getLower(3); l<row.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", row.getGenIndex(idx) == count++);
        }
}

void
test_5()
{
  int start[2] = {1, 2};
  unsigned shape[2] = {10, 15};
  Lucee::RowMajorIndexer<2> row(shape, start);
  Lucee::RowMajorIndexer<2> rowCopy(row); // copy it
 
// now test the copied indexer (testing copy ctor)
  int count=0;
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
      LC_ASSERT("Testing 2D indexer", rowCopy.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", rowCopy.getGenIndex(idx) == count++);
    }
}

void
test_6()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::RowMajorIndexer<4> row(shape, start);
  Lucee::RowMajorIndexer<4> rowCopy(row);

  int count=0;
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
      for (int k=rowCopy.getLower(2); k<rowCopy.getUpper(2); ++k)
        for (int l=rowCopy.getLower(3); l<rowCopy.getUpper(3); ++l)
          LC_ASSERT("Testing 4D indexer", rowCopy.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
      for (int k=rowCopy.getLower(2); k<rowCopy.getUpper(2); ++k)
        for (int l=rowCopy.getLower(3); l<rowCopy.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", rowCopy.getGenIndex(idx) == count++);
        }
}

void
test_7()
{
  int start[4] = {1, 2, -2, 3};
  unsigned shape[4] = {10, 15, 20, 12};
  Lucee::RowMajorIndexer<4> row(shape, start);

  int zeros[4] = {0, 0, 0, 0};
  unsigned ones[4] = {1, 1, 1, 1};
  Lucee::RowMajorIndexer<4> rowCopy(ones, zeros);
  rowCopy = row;

  int count=0;
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
      for (int k=rowCopy.getLower(2); k<rowCopy.getUpper(2); ++k)
        for (int l=rowCopy.getLower(3); l<rowCopy.getUpper(3); ++l)
          LC_ASSERT("Testing 4D indexer", rowCopy.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int i=rowCopy.getLower(0); i<rowCopy.getUpper(0); ++i)
    for (int j=rowCopy.getLower(1); j<rowCopy.getUpper(1); ++j)
      for (int k=rowCopy.getLower(2); k<rowCopy.getUpper(2); ++k)
        for (int l=rowCopy.getLower(3); l<rowCopy.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", rowCopy.getGenIndex(idx) == count++);
        }
}

void
test_1_1()
{
  int lower[1] = {1};
  int upper[1] = {10+1};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<1> row(rgn);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    LC_ASSERT("Testing 1D indexer", row.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indexer", row.getGenIndex(idx) == count++);
  }
}

void
test_2_2()
{
  int lower[2] = {1, 2};
  int upper[2] = {10+1, 15+2};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<2> row(rgn);

  int count=0;
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      LC_ASSERT("Testing 2D indexer", row.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indexer", row.getGenIndex(idx) == count++);
    }
}

void
test_3_3()
{
  int lower[3] = {1, 2, -2};
  int upper[3] = {10+1, 15+2, 20-2};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<3> row(rgn);

  int count=0;

  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        LC_ASSERT("Testing 3D indexer", row.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indexer", row.getGenIndex(idx) == count++);
      }
}

void
test_4_4()
{
  int lower[4] = {1, 2, -2, 3};
  int upper[4] = {10+1, 15+2, 20-2, 12+3};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<4> row(rgn);

  int count=0;

  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        for (int l=row.getLower(3); l<row.getUpper(3); ++l)
          LC_ASSERT("Testing 4D indexer", row.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int i=row.getLower(0); i<row.getUpper(0); ++i)
    for (int j=row.getLower(1); j<row.getUpper(1); ++j)
      for (int k=row.getLower(2); k<row.getUpper(2); ++k)
        for (int l=row.getLower(3); l<row.getUpper(3); ++l)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indexer", row.getGenIndex(idx) == count++);
        }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcrowmajorindexer");
  test_1();
  test_2();
  test_3();
  test_4();
  test_5();
  test_6();
  test_7();

  test_1_1();
  test_2_2();
  test_3_3();
  test_4_4();
  LC_END_TESTS;
}
