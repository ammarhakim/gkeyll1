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
test_array_1()
{
  unsigned shape[2] = {5, 10};
  Lucee::Array<2, double> arr1(shape);

  int start[2] = {2, 0};
  Lucee::Array<2, double> arr2(shape, start);

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking start indices", arr1.getStart(i) == 0);
    LC_ASSERT("Checking start indices", arr2.getStart(i) == start[i]);

    LC_ASSERT("Checking end indices", arr1.getEnd(i) == (int) shape[i]);
    int end = start[i]+shape[i];
    LC_ASSERT("Checking end indices", arr2.getEnd(i) == end);

    LC_ASSERT("Checking shape", arr1.getShape(i) == shape[i]);
  }
    
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcarray");
  test_array_1();
  LC_END_TESTS;
}
