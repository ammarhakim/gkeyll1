/**
 * @file	lcfield.cxx
 *
 * @brief	Unit tests for Lucee::Field class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcField.h>
#include <LcTest.h>

void
test_1()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcField(rgn, 3, 10.0);

  LC_ASSERT("Testing if number of components is correct", 
    elcField.getNumComponents() == 3);

  Lucee::Region<2, int> idxRgn = elcField.getRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      idxRgn.getLower(i) == rgn.getLower(i));
    LC_ASSERT("Testing if index region is correct",
      idxRgn.getUpper(i) == rgn.getUpper(i));
  }

  for (int i=elcField.getLower(0); i<elcField.getUpper(0); ++i)
    for (int j=elcField.getLower(1); j<elcField.getUpper(1); ++j)
    {
      for (int k=0; k<elcField.getNumComponents(); ++k)
        LC_ASSERT("Testing default values", elcField(i,j,k) == 10.0);
    }

  elcField = 3.0;
  for (int i=elcField.getLower(0); i<elcField.getUpper(0); ++i)
    for (int j=elcField.getLower(1); j<elcField.getUpper(1); ++j)
    {
      for (int k=0; k<elcField.getNumComponents(); ++k)
        LC_ASSERT("Testing assigned values", elcField(i,j,k) == 3.0);
    }
}

void
test_2()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcField(rgn, 3, 10.0);

// create ptr to field elements
  Lucee::FieldPtr<double> ptr = elcField.createPtr();

  LC_ASSERT("Testing if numComponents is 3", ptr.getNumComponents() == 3);

// set values at initial location (0,0)
  ptr[0] = 1.5;
  ptr[1] = 2.5;
  ptr[2] = 3.5;

// see if setting worked
  LC_ASSERT("Testing if assignment worked", ptr[0] == 1.5);
  LC_ASSERT("Testing if assignment worked", ptr[1] == 2.5);
  LC_ASSERT("Testing if assignment worked", ptr[2] == 3.5);

  LC_ASSERT("Testing if assignment worked", elcField(0,0,0) == 1.5);
  LC_ASSERT("Testing if assignment worked", elcField(0,0,1) == 2.5);
  LC_ASSERT("Testing if assignment worked", elcField(0,0,2) == 3.5);

  elcField.setPtr(ptr, 5, 5);
  ptr[0] = 1.5;
  ptr[1] = 2.5;
  ptr[2] = 3.5;

// see if setting worked
  LC_ASSERT("Testing if assignment worked", ptr[0] == 1.5);
  LC_ASSERT("Testing if assignment worked", ptr[1] == 2.5);
  LC_ASSERT("Testing if assignment worked", ptr[2] == 3.5);

  LC_ASSERT("Testing if assignment worked", elcField(5,5,0) == 1.5);
  LC_ASSERT("Testing if assignment worked", elcField(5,5,1) == 2.5);
  LC_ASSERT("Testing if assignment worked", elcField(5,5,2) == 3.5);
}

void
test_3()
{
  int lower[1] = {2};
  int upper[1] = {12};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::Field<1, double> elcFld(rgn, 3, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
  {
    elcFld.setPtr(ptr, i);
    for (unsigned n=0; n<ptr.getNumComponents(); ++n)
      ptr[n] = (i+0.5)*n;
  }

  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
  {
    for (int k=0; k<elcFld.getNumComponents(); ++k)
      LC_ASSERT("Testing if 1D setting worked", elcFld(i,k) == (i+0.5)*k);
  }
}

void
test_4()
{
  int lower[2] = {2, 5};
  int upper[2] = {12, 32};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcFld(rgn, 2, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<rgn.getUpper(1); ++j)
    {
      elcFld.setPtr(ptr, i, j);
      for (unsigned n=0; n<ptr.getNumComponents(); ++n)
        ptr[n] = (i+3*j+0.5)*n;
    }

  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<rgn.getUpper(1); ++j)
    {
      for (int k=0; k<elcFld.getNumComponents(); ++k)
        LC_ASSERT("Testing if 2D setting worked", elcFld(i,j,k) == (i+3*j+0.5)*k);
    }
}

void
test_5()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> elcFld(rgn, 3, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<rgn.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<rgn.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

  for (int i=elcFld.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<rgn.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<rgn.getUpper(2); ++k)
      {
        for (int n=0; n<elcFld.getNumComponents(); ++n)
          LC_ASSERT("Testing if 3D setting worked", elcFld(i,j,k,n) == (i+3*j+5.5*k+0.5)*n);
      }
}

int
main(void)
{
  LC_BEGIN_TESTS("lcfield");
  test_1();
  test_2();
  test_3();
  test_4();
  test_5();
  LC_END_TESTS;
}
