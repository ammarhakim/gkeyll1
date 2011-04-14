/**
 * @file	lcvec3.cxx
 *
 * @brief	Unit tests for Lucee::FixedVector class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcTest.h>
#include <LcVec3.h>

// std includes
#include <cmath>

void
test_1()
{
  Lucee::Vec3<double> x(1,0,0), y(0,1,0), z(0,0,1);

  LC_ASSERT("Testing dot products", x.dot(x) == 1);
  LC_ASSERT("Testing dot products", x.dot(y) == 0);
  LC_ASSERT("Testing dot products", x.dot(z) == 0);
  LC_ASSERT("Testing dot products", y.dot(y) == 1);
  LC_ASSERT("Testing dot products", y.dot(z) == 0);
  LC_ASSERT("Testing dot products", z.dot(z) == 1);

  Lucee::Vec3<double> xy = x+y;
  LC_ASSERT("Testing sum of vectors", xy[0] == 1);
  LC_ASSERT("Testing sum of vectors", xy[1] == 1);
  LC_ASSERT("Testing sum of vectors", xy[2] == 0);

  Lucee::Vec3<double> xyz = x+y+z;
  LC_ASSERT("Testing sum of vectors", xyz[0] == 1);
  LC_ASSERT("Testing sum of vectors", xyz[1] == 1);
  LC_ASSERT("Testing sum of vectors", xyz[2] == 1);

  Lucee::Vec3<double> xy1 = x-y;
  LC_ASSERT("Testing sum of vectors", xy1[0] == 1);
  LC_ASSERT("Testing sum of vectors", xy1[1] == -1);
  LC_ASSERT("Testing sum of vectors", xy1[2] == 0);

  Lucee::Vec3<double> xyz1 = x-y-z;
  LC_ASSERT("Testing sum of vectors", xyz1[0] == 1);
  LC_ASSERT("Testing sum of vectors", xyz1[1] == -1);
  LC_ASSERT("Testing sum of vectors", xyz1[2] == -1);

  xyz.normalize();
  LC_ASSERT("Testing if normalization worked", xyz.getNorm() == 1.0);
  LC_ASSERT("Testing if normalization worked", xyz[0] == 1/std::sqrt(3));
  LC_ASSERT("Testing if normalization worked", xyz[1] == 1/std::sqrt(3));
  LC_ASSERT("Testing if normalization worked", xyz[2] == 1/std::sqrt(3));

  Lucee::Vec3<double> xy3 = x.cross(y);
  LC_ASSERT("Testing if cross product", xy3[0] == 0);
  LC_ASSERT("Testing if cross product", xy3[1] == 0);
  LC_ASSERT("Testing if cross product", xy3[2] == 1);

  Lucee::Vec3<double> yz3 = y.cross(z);
  LC_ASSERT("Testing if cross product", yz3[0] == 1);
  LC_ASSERT("Testing if cross product", yz3[1] == 0);
  LC_ASSERT("Testing if cross product", yz3[2] == 0);

  Lucee::Vec3<double> zx3 = z.cross(x);
  LC_ASSERT("Testing if cross product", zx3[0] == 0);
  LC_ASSERT("Testing if cross product", zx3[1] == 1);
  LC_ASSERT("Testing if cross product", zx3[2] == 0);

  xyz = x+y+z;
  xyz.scale(0.5);
  LC_ASSERT("Testing sum of vectors", xyz[0] == 0.5);
  LC_ASSERT("Testing sum of vectors", xyz[1] == 0.5);
  LC_ASSERT("Testing sum of vectors", xyz[2] == 0.5);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcvec3");
  test_1();
  LC_END_TESTS;
}
