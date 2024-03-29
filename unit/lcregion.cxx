/**
 * @file	lcregion.cxx
 *
 * @brief	Unit tests for Lucee::Region class
 */

// lucee includes
#include <LcRegion.h>
#include <LcRowMajorSequencer.h>
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
  LC_ASSERT("Testing containment", ibox.contains(ibox1) == true);

  int lo2[2] = {3, 4};
  int up2[2] = {14, 10};
  Lucee::Region<2, int> ibox2(lo2, up2);
  LC_ASSERT("Testing containment", ibox.contains(ibox2) == false);

  int lo3[2] = {0, 0};
  int up3[2] = {10, 10};
  Lucee::Region<2, int> ibox3(lo3, up3);

  int lo4[2] = {5, 5};
  int up4[2] = {15, 15};
  Lucee::Region<2, int> ibox4(lo4, up4);

  Lucee::Region<2, int> ibox34 = ibox3.intersect(ibox4);

  LC_ASSERT("Checking intersected box", ibox34.getLower(0) == 5);
  LC_ASSERT("Checking intersected box", ibox34.getLower(1) == 5);

  LC_ASSERT("Checking intersected box", ibox34.getUpper(0) == 10);
  LC_ASSERT("Checking intersected box", ibox34.getUpper(1) == 10);

  LC_ASSERT("Checking volume of intersection", ibox34.getVolume() == 25);
}

void
test_9()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox(lower, upper);

  Lucee::Region<2, int> ibox0 = ibox.deflate(0);
  Lucee::Region<2, int> ibox1 = ibox.deflate(1);
  LC_ASSERT("Testing volume of deflated boxes", ibox0.getVolume() == 10);
  LC_ASSERT("Testing volume of deflated boxes", ibox1.getVolume() == 10);

  Lucee::RowMajorSequencer<2> seq0(ibox0);
  Lucee::RowMajorSequencer<2> seq1(ibox1);

  int j=ibox.getLower(1);
  while (seq0.step())
  {
    int idx[2];
    seq0.fillWithIndex(idx);
    for (int i=ibox.getLower(0); i<ibox.getUpper(0); ++i)
    {
      idx[0] = i;
      LC_ASSERT("Testing deflated box indices", idx[0]==i && idx[1]==j);
    }
    j++;
  }
  
  int i=ibox.getLower(0);
  while (seq1.step())
  {
    int idx[2];
    seq1.fillWithIndex(idx);
    for (int j=ibox.getLower(1); j<ibox.getUpper(1); ++j)
    {
      idx[1] = j;
      LC_ASSERT("Testing deflated box indices", idx[0]==i && idx[1]==j);
    }
    i++;
  }
}

void
test_10()
{  
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox(lower, upper);

  LC_ASSERT("Checking if box is empty", ibox.isEmpty() == false);

  Lucee::Region<2, int> ebox;
  LC_ASSERT("Checking if box is empty", ebox.isEmpty() == true);
}

void
test_11()
{  
  int lower[2] = {3, 4};
  int upper[2] = {13, 14};
  Lucee::Region<2, int> ibox(lower, upper), ebox(lower, upper);
  upper[1] = 15;
  Lucee::Region<2, int> nebox(lower, upper);

  LC_ASSERT("Compare boxes", ibox == ebox);
  LC_ASSERT("Compare boxes", (ibox == nebox) == false);

  LC_ASSERT("Compare boxes", ibox != nebox);
  LC_ASSERT("Compare boxes", (ibox != ebox) == false);
}

void
test_12()
{
  int lower[2] = {3, 4};
  int upper[2] = {13, 20};
  Lucee::Region<2, int> ibox(lower, upper);

  Lucee::Region<2, int> ibox0 = ibox.resetBounds(0, 1, 4);
  LC_ASSERT("Testing resetBound method", ibox0.getLower(0) == 1);
  LC_ASSERT("Testing resetBound method", ibox0.getUpper(0) == 4);

  Lucee::Region<2, int> ibox1 = ibox.resetBounds(1, 15, 17);
  LC_ASSERT("Testing resetBound method", ibox1.getLower(1) == 15);
  LC_ASSERT("Testing resetBound method", ibox1.getUpper(1) == 17);
}

void
test_13()
{
  int lower[2] = {1, 0};
  int upper[2] = {10, 12};
  Lucee::Region<2, int> rgn(lower, upper);

  LC_ASSERT("Testing volume before reset", rgn.getVolume() == 9*12);

  rgn.setLower(0, 3);
  LC_ASSERT("Testing new bounds", rgn.getLower(0) == 3);
  LC_ASSERT("Testing volume after reset", rgn.getVolume() == (10-3)*(12-0));

  rgn.setUpper(1, 15);
  LC_ASSERT("Testing new bounds", rgn.getUpper(1) == 15);
  LC_ASSERT("Testing volume after reset", rgn.getVolume() == (10-3)*(15-0));
}

int
main(int argc, char **argv) 
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
  test_9();
  test_10();
  test_11();
  test_12();
  test_13();
  LC_END_TESTS;
}
