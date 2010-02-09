/**
 * @file	lcregion.cxx
 *
 * @brief	Unit tests for Lucee::Region class
 *
 * @version	$Id: lcregion.cxx 162 2009-08-28 20:07:10Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcRegion.h>
#include <LcTest.h>

void
test_1()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox(lower, upper);
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking lower bounds", ibox.getLower(i) == lower[i]);
    LC_ASSERT("Checking upper bounds", ibox.getUpper(i) == upper[i]);
    LC_ASSERT("Checking shape", ibox.getShape(i) == (upper[i]-lower[i]));
  }
  LC_ASSERT("Checking volume", ibox.getVolume() == 100);

  Lucee::FixedVector<2, int> lfv = ibox.getLower();
  Lucee::FixedVector<2, int> ufv = ibox.getUpper();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking lower bounds", ibox.getLower(i) == lfv[i]);
    LC_ASSERT("Checking upper bounds", ibox.getUpper(i) == ufv[i]);
  }

  int shape[2] = {10, 12};
  Lucee::Region<2, int> ibox2(shape);
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking lower bounds", ibox2.getLower(i) == 0);
    LC_ASSERT("Checking upper bounds", ibox2.getUpper(i) == shape[i]);
    LC_ASSERT("Checking shape", ibox2.getShape(i) == shape[i]);
  }
  LC_ASSERT("Checking volume", ibox2.getVolume() == 120);

  Lucee::FixedVector<2, int> lfv2 = ibox2.getLower();
  Lucee::FixedVector<2, int> ufv2 = ibox2.getUpper();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking lower bounds", ibox2.getLower(i) == lfv2[i]);
    LC_ASSERT("Checking upper bounds", ibox2.getUpper(i) == ufv2[i]);
  }
}

void
test_2()
{
  double lower[2] = {0, 0};
  double upper[2] = {1.0, 1.0};
  Lucee::Region<2, double> rgn(lower, upper);

  double p1[2] = {0.5, 0.5};
  LC_ASSERT("Checking if point is inside region",
    rgn.isInside(p1) == true);

  double p2[2] = {0.0, 0.5};
  LC_ASSERT("Checking if point is inside region",
    rgn.isInside(p2) == true);

  double p3[2] = {-1.0, 0.5};
  LC_ASSERT("Checking if point is inside region",
    rgn.isInside(p3) == false);

  double p4[2] = {0.5, 1.0};
  LC_ASSERT("Checking if point is inside region",
    rgn.isInside(p4) == false);

  double p5[2] = {1.0, 0.5};
  LC_ASSERT("Checking if point is inside region",
    rgn.isInside(p5) == false);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcregion");
  test_1();
  test_2();
  LC_END_TESTS;
}
