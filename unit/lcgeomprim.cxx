/**
 * @file	lcgeomprim.cxx
 *
 * @brief	Unit tests for geometry primitives.
 */

// lucee includes
#include <LcMathLib.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  Lucee::Vec3<double> a(0, 0, 0), b(1, 0, 0), c(1, 1, 0);
  LC_ASSERT("Testing area of triangle", Lucee::calcTriArea(a, b, c) == 0.5);

  c[1] = 0.5; // c = (1, 0.5, 0)
  LC_ASSERT("Testing area of triangle", Lucee::calcTriArea(a, b, c) == 0.25);

  c[0] = 0.5;
  c[1] = 0.5*std::sqrt(3.0); // c = (0.5, sqrt(3)/2, 0)
  LC_ASSERT("Testing area of triangle", Lucee::calcTriArea(a, b, c) == 0.5*0.5*std::sqrt(3.0));
}

void
test_2()
{
  Lucee::Vec3<double> a(0, 0, 0), b(1, 0, 0), c(1, 1, 0), d(0, 1, 0);
  LC_ASSERT("Testing area of quad", Lucee::calcQuadArea(a, b, c, d) == 1.0);
}

int
main (int argc, char **argv)
{
  LC_BEGIN_TESTS("lcgeomprim");
  test_1();
  test_2();
  LC_END_TESTS;
}

