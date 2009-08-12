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
    LC_ASSERT("Testing 1D indedxer", col.getIndex(i) == count++);

  count=0;
  int idx[1];
  for (int i=col.getLower(0); i<col.getUpper(0); ++i)
  {
    idx[0] = i;
    LC_ASSERT("Testing 1D indedxer", col.getGenIndex(idx) == count++);
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
      LC_ASSERT("Testing 2D indedxer", col.getIndex(i,j) == count++);

  count=0;
  int idx[2];
  for (int j=col.getLower(1); j<col.getUpper(1); ++j)
    for (int i=col.getLower(0); i<col.getUpper(0); ++i)
    {
      idx[0] = i; idx[1] = j;
      LC_ASSERT("Testing 2D indedxer", col.getGenIndex(idx) == count++);
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
        LC_ASSERT("Testing 3D indedxer", col.getIndex(i,j,k) == count++);

  count=0;
  int idx[3];
  for (int k=col.getLower(2); k<col.getUpper(2); ++k)
    for (int j=col.getLower(1); j<col.getUpper(1); ++j)
      for (int i=col.getLower(0); i<col.getUpper(0); ++i)
      {
        idx[0] = i; idx[1] = j; idx[2] = k;
        LC_ASSERT("Testing 3D indedxer", col.getGenIndex(idx) == count++);
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
          LC_ASSERT("Testing 4D indedxer", col.getIndex(i,j,k,l) == count++);

  count=0;
  int idx[4];
  for (int l=col.getLower(3); l<col.getUpper(3); ++l)
    for (int k=col.getLower(2); k<col.getUpper(2); ++k)
      for (int j=col.getLower(1); j<col.getUpper(1); ++j)
        for (int i=col.getLower(0); i<col.getUpper(0); ++i)
        {
          idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
          LC_ASSERT("Testing 4D indedxer", col.getGenIndex(idx) == count++);
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
  LC_END_TESTS;
}
