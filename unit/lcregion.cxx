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

void
test_3()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox_o(lower, upper);

  Lucee::Region<2, int> ibox(ibox_o);
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
}

void
test_4()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox_o(lower, upper);

  Lucee::Region<2, int> ibox(lower, upper);
  ibox = ibox_o;
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
}

void
test_5()
{
  int lower[3] = {1, 1, 1};
  int upper[3] = {10, 10, 10};
  Lucee::Region<3, int> rgn(lower, upper);

  int lowerExt[3] = {1, 1, 1};
  int upperExt[3] = {3, 3, 3};
  Lucee::Region<3, int> extRgn(rgn.extend(lowerExt, upperExt));
  
  for (unsigned i=0; i<3; ++i)
  {
    LC_ASSERT("Checking extended region", extRgn.getLower(i) == 0);
    LC_ASSERT("Checking extended region", extRgn.getUpper(i) == 13);
  }
  LC_ASSERT("Checking extended region volume", extRgn.getVolume() == 13*13*13);
}

void
test_6()
{
  int lower[2] = {1, 0};
  int upper[2] = {10, 12};
  Lucee::Region<2, int> rgn(lower, upper);

  Lucee::Region<3, int> infRgn(rgn.inflate(-2, 10));
  
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Checking inflated region lower bounds", 
      infRgn.getLower(i) == rgn.getLower(i));
    LC_ASSERT("Checking inflated region upper bounds", 
      infRgn.getUpper(i) == rgn.getUpper(i));
  }
  LC_ASSERT("Checking inflated region lower bounds", 
    infRgn.getLower(2) == -2);
  LC_ASSERT("Checking inflated region upper bounds", 
    infRgn.getUpper(2) == 10);

  LC_ASSERT("Checking inflated region volume", 
    infRgn.getVolume() == rgn.getVolume()*(10+2));
}

void
test_7()
{
  int lower[2] = {1,2};
  int upper[2] = {10, 10};
  Lucee::Region<2, int> rgn(lower, upper);

  int lowerExp[2] = {1, 2};
  int upperExp[2] = {1, 2};
  Lucee::Region<2, int> expRgn = rgn.extend(lowerExp, upperExp);

  LC_ASSERT("Checking extended box", expRgn.getLower(0) == 0);
  LC_ASSERT("Checking extended box", expRgn.getLower(1) == 0);
  LC_ASSERT("Checking extended box", expRgn.getUpper(0) == 11);
  LC_ASSERT("Checking extended box", expRgn.getUpper(1) == 12);
}

void
test_8()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox(lower, upper);

  LC_ASSERT("Testing self containment", ibox.contains(ibox) == true);

  int lo1[2] = {3, 4};
  int up1[2] = {10, 10};
  Lucee::Region<2, int> ibox1(lo1, up1);
  //LC_ASSERT("Testing self containment", ibox.contains(ibox1) == true);

  int lo2[2] = {3, 4};
  int up2[2] = {14, 10};
  Lucee::Region<2, int> ibox2(lo2, up2);
  LC_ASSERT("Testing self containment", ibox.contains(ibox2) == false);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcregion");
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
