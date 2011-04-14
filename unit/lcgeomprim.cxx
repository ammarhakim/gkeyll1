/**
 * @file	lcgeomprim.cxx
 *
 * @brief	Unit tests for geometry primitives.
 *
 * @version	$Id$
 */

// lucee includes
#include <LcMathLib.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::Vec3<double> a(0, 0, 0), b(0, 1, 0), c(1, 0, 0);
  LC_ASSERT("Testing area of triangle", Lucee::calcTriArea(a, b, c) == 0.5);
  LC_ASSERT("Testing area of triangle", Lucee::calcTriArea(b, a, c) == 0.5);
}

int
main (void)
{
  LC_BEGIN_TESTS("lcgeomprim");
  test_1();
  LC_END_TESTS;
}

