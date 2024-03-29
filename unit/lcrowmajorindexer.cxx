/**
 * @file	lcrowmajorindexer.cxx
 *
 * @brief	Unit tests for Lucee::RowMajorIndexer class
 */

// lucee includes
#include <LcRowMajorIndexer.h>
#include <LcRowMajorSequencer.h>
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
    LC_ASSERT("Testing 1D indexer", row.getIndex(idx) == count++);
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
      LC_ASSERT("Testing 2D indexer", row.getIndex(idx) == count++);
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
        LC_ASSERT("Testing 3D indexer", row.getIndex(idx) == count++);
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
          LC_ASSERT("Testing 4D indexer", row.getIndex(idx) == count++);
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
      LC_ASSERT("Testing 2D indexer", rowCopy.getIndex(idx) == count++);
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
          LC_ASSERT("Testing 4D indexer", rowCopy.getIndex(idx) == count++);
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
          LC_ASSERT("Testing 4D indexer", rowCopy.getIndex(idx) == count++);
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
    LC_ASSERT("Testing 1D indexer", row.getIndex(idx) == count++);
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
      LC_ASSERT("Testing 2D indexer", row.getIndex(idx) == count++);
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
        LC_ASSERT("Testing 3D indexer", row.getIndex(idx) == count++);
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
          LC_ASSERT("Testing 4D indexer", row.getIndex(idx) == count++);
        }
}

void
test_8()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<4> row(rgn);

// deflate the 4D indexer into a 2D indexer
  unsigned defDims[2] = {1, 2};
  int defDimsIdx[2] = {10, 10};

  Lucee::RowMajorIndexer<2> defRow
    = row.deflate<2>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defRow.getLower(0) == 3);
  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defRow.getLower(1) == 10);

  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defRow.getUpper(0) == 15);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defRow.getUpper(1) == 30);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      LC_ASSERT("Testing deflated indexer", defRow.getIndex(i,j) == row.getIndex(i,10,10,j));
}

void
test_9()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<4> row(rgn);

// deflate the 4D indexer into a 2D indexer
  unsigned defDims[2] = {0, 3};
  int defDimsIdx[2] = {4, 15};

  Lucee::RowMajorIndexer<2> defRow
    = row.deflate<2>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defRow.getLower(0) == 6);
  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defRow.getLower(1) == 8);

  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defRow.getUpper(0) == 20);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defRow.getUpper(1) == 25);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      LC_ASSERT("Testing deflated indexer", defRow.getIndex(i,j) == row.getIndex(4,i,j,15));
}

void
test_10()
{
  int lower[4] = {3, 6, 8, 10};
  int upper[4] = {15, 20, 25, 30};
  Lucee::Region<4, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<4> row(rgn);

// deflate the 4D indexer into a 1D indexer
  unsigned defDims[3] = {0, 1, 2};
  int defDimsIdx[3] = {4, 15, 20};

  Lucee::RowMajorIndexer<1> defRow
    = row.deflate<1>(defDims, defDimsIdx);

  LC_ASSERT("Testing lower bounds of deflated indexer", 
    defRow.getLower(0) == 10);
  LC_ASSERT("Testing upper bounds of deflated indexer", 
    defRow.getUpper(0) == 30);

// now check deflated indexer
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    LC_ASSERT("Testing deflated indexer", defRow.getIndex(i) == row.getIndex(4,15,20,i));
}

void
test_11()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 20};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::RowMajorIndexer<2> rowIdxr(rgn);

  int finalIdx = rowIdxr.getIndex(rgn.getUpper(0)-1, rgn.getUpper(1)-1);
  LC_ASSERT("Testing indexing beyond bounds", rowIdxr.getIndex(rgn.getUpper(0), rgn.getUpper(1)-1) > finalIdx);
  LC_ASSERT("Testing indexing beyond bounds", rowIdxr.getIndex(rgn.getLower(0)-1, rgn.getLower(1)) < 0);
}

void
test_12()
{
  int zeros[3], ishape[3];
  unsigned shape[3];
  for (unsigned d=0; d<3; ++d)
    zeros[d] = 0;
  shape[0] = 10; shape[1] = 3; shape[2] = 7;
  for (unsigned d=0; d<3; ++d)
    ishape[d] = shape[d]; // I SOMETIMES HATE C++

  Lucee::RowMajorIndexer<3> cIdx(shape, zeros);
  Lucee::Region<3, int> rgn(ishape);
  Lucee::RowMajorSequencer<3> cSeq(rgn);

  int idx[3];
  unsigned count = 0;
  while (cSeq.step())
  {
    cSeq.fillWithIndex(idx);
    LC_ASSERT("Checking if indexer/sequencer match", cIdx.getIndex(idx) == count);
    count++;
  }
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcrowmajorindexer");
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

  test_1_1();
  test_2_2();
  test_3_3();
  test_4_4();
  LC_END_TESTS;
}
