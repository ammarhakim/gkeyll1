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
#include <LcVector.h>

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
      for (unsigned k=0; k<elcField.getNumComponents(); ++k)
        LC_ASSERT("Testing default values", elcField(i,j,k) == 10.0);
    }

  elcField = 3.0;
  for (int i=elcField.getLower(0); i<elcField.getUpper(0); ++i)
    for (int j=elcField.getLower(1); j<elcField.getUpper(1); ++j)
    {
      for (unsigned k=0; k<elcField.getNumComponents(); ++k)
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
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
  {
    elcFld.setPtr(ptr, i);
    for (unsigned n=0; n<ptr.getNumComponents(); ++n)
      ptr[n] = (i+0.5)*n;
  }

  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
  {
    for (unsigned k=0; k<elcFld.getNumComponents(); ++k)
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
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
    {
      elcFld.setPtr(ptr, i, j);
      for (unsigned n=0; n<ptr.getNumComponents(); ++n)
        ptr[n] = (i+3*j+0.5)*n;
    }

  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
    {
      for (unsigned k=0; k<elcFld.getNumComponents(); ++k)
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
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        for (unsigned n=0; n<elcFld.getNumComponents(); ++n)
          LC_ASSERT("Testing if 3D setting worked", elcFld(i,j,k,n) == (i+3*j+5.5*k+0.5)*n);
      }
}

void
test_6()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> elcFld(rgn, 3, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

// create constant reference to field
  const Lucee::Field<3, double>& elcFldCnst = elcFld;
  Lucee::ConstFieldPtr<double> cPtr = elcFldCnst.createConstPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(cPtr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          LC_ASSERT("Testing if setting worked", cPtr[n] == (i+3*j+5.5*k+0.5)*n);
      }
}

void
test_7()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> elcFld(rgn, 3, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

  int vlo[3] = {4, 8, 8};
  int vup[3] = {10, 30, 10};
  Lucee::Region<3, int> vrgn(vlo, vup);
  Lucee::Field<3, double> subElcFld = elcFld.getView(vrgn);

// check bounds of view
  for (unsigned i=0; i<3; ++i)
  {
    LC_ASSERT("Testing bounds of view", subElcFld.getLower(i) == vlo[i]);
    LC_ASSERT("Testing bounds of view", subElcFld.getUpper(i) == vup[i]);
  }
  LC_ASSERT("Testing number of components", subElcFld.getNumComponents() == 3);

  for (int i=subElcFld.getLower(0); i<subElcFld.getUpper(0); ++i)
    for (int j=subElcFld.getLower(1); j<subElcFld.getUpper(1); ++j)
      for (int k=subElcFld.getLower(2); k<subElcFld.getUpper(2); ++k)
        for (int n=0; n<3; ++n)
          LC_ASSERT("Testing view components", subElcFld(i,j,k,n) == (i+3*j+5.5*k+0.5)*n);

  Lucee::FieldPtr<double> subPtr = subElcFld.createPtr();
  for (int i=subElcFld.getLower(0); i<subElcFld.getUpper(0); ++i)
    for (int j=subElcFld.getLower(1); j<subElcFld.getUpper(1); ++j)
      for (int k=subElcFld.getLower(2); k<subElcFld.getUpper(2); ++k)
      {
        subElcFld.setPtr(subPtr, i, j, k);
        for (unsigned n=0; n<subPtr.getNumComponents(); ++n)
          LC_ASSERT("Testing view components", subPtr[n] == (i+3*j+5.5*k+0.5)*n);
      }
}

void
test_8()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> emFld(rgn, 6, 10.0);

  Lucee::FieldPtr<double> ptr = emFld.createPtr();
  for (int i=emFld.getLower(0); i<emFld.getUpper(0); ++i)
    for (int j=emFld.getLower(1); j<emFld.getUpper(1); ++j)
      for (int k=emFld.getLower(2); k<emFld.getUpper(2); ++k)
      {
        emFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

// now extract first three components of field
  Lucee::Field<3, double> elcFld = emFld.getSubCompView(0, 3);

  LC_ASSERT("Testing number of components in sub-comp view", elcFld.getNumComponents() == 3);
  Lucee::FieldPtr<double> elcPtr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(elcPtr, i, j, k);
        for (unsigned n=0; n<elcPtr.getNumComponents(); ++n)
          LC_ASSERT("Testing components of sub-comp view", elcPtr[n] == (i+3*j+5.5*k+0.5)*n);
      }

// now extract last three components of field
  Lucee::Field<3, double> magFld = emFld.getSubCompView(3, 6);

  LC_ASSERT("Testing number of components in sub-comp view", magFld.getNumComponents() == 3);
  LC_ASSERT("Testing components bounds", magFld.getLower(3) == 0);
  LC_ASSERT("Testing components bounds", magFld.getUpper(3) == 3);

  Lucee::FieldPtr<double> magPtr = magFld.createPtr();
  for (int i=magFld.getLower(0); i<magFld.getUpper(0); ++i)
    for (int j=magFld.getLower(1); j<magFld.getUpper(1); ++j)
      for (int k=magFld.getLower(2); k<magFld.getUpper(2); ++k)
      {
        magFld.setPtr(magPtr, i, j, k);
        for (unsigned n=0; n<magPtr.getNumComponents(); ++n)
          LC_ASSERT("Testing components of sub-comp view", magPtr[n] == (i+3*j+5.5*k+0.5)*(n+3));
      }

  for (int i=magFld.getLower(0); i<magFld.getUpper(0); ++i)
    for (int j=magFld.getLower(1); j<magFld.getUpper(1); ++j)
      for (int k=magFld.getLower(2); k<magFld.getUpper(2); ++k)
        for (int n=magFld.getLower(3); n<magFld.getUpper(3); ++n)
          LC_ASSERT("Testing components of sub-comp view", 
            magFld(i,j,k,n) == (i+3*j+5.5*k+0.5)*(n+3));
}

void
test_9()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> elcFld(rgn, 3, 10.0);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }

  int vlo[3] = {4, 8, 8};
  int vup[3] = {10, 30, 10};
  Lucee::Region<3, int> vrgn(vlo, vup);
  Lucee::Field<3, double> subElcFld = elcFld.getView(vrgn);

// now select last two components of view-field
  Lucee::Field<3, double> eyez = subElcFld.getSubCompView(1,3);

  LC_ASSERT("Testing numComponents of sub-comp field", eyez.getNumComponents() == 2);

  Lucee::ConstFieldPtr<double> eyezPtr = eyez.createConstPtr();
  for (int i=eyez.getLower(0); i<eyez.getUpper(0); ++i)
    for (int j=eyez.getLower(1); j<eyez.getUpper(1); ++j)
      for (int k=eyez.getLower(2); k<eyez.getUpper(2); ++k)
      {
        eyez.setPtr(eyezPtr, i, j, k);
        for (unsigned n=0; n<eyez.getNumComponents(); ++n)
          LC_ASSERT("Testing sub-comp field", eyezPtr[n] == (i+3*j+5.5*k+0.5)*(n+1));
      }
}

void
test_10()
{
  int lower[3] = {2, 5, 6};
  int upper[3] = {12, 32, 12};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::Field<3, double> elcFld(rgn, 3);

  Lucee::FieldPtr<double> ptr = elcFld.createPtr();
  for (int i=elcFld.getLower(0); i<elcFld.getUpper(0); ++i)
    for (int j=elcFld.getLower(1); j<elcFld.getUpper(1); ++j)
      for (int k=elcFld.getLower(2); k<elcFld.getUpper(2); ++k)
      {
        elcFld.setPtr(ptr, i, j, k);
        for (unsigned n=0; n<ptr.getNumComponents(); ++n)
          ptr[n] = (i+3*j+5.5*k+0.5)*n;
      }
}

void
test_11()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  int lg[2] = {1, 3};
  int ug[2] = {2, 5};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcField(rgn, 3, lg, ug, 10.0);

  LC_ASSERT("Testing if number of components is correct", 
    elcField.getNumComponents() == 3);

  Lucee::Region<2, int> localRgn = elcField.getRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      localRgn.getLower(i) == rgn.getLower(i));
    LC_ASSERT("Testing if index region is correct",
      localRgn.getUpper(i) == rgn.getUpper(i));
  }

  Lucee::Region<2, int> extRgn = elcField.getExtRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      extRgn.getLower(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      extRgn.getUpper(i) == rgn.getUpper(i)+ug[i]);
  }

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      elcField.getLowerExt(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      elcField.getUpperExt(i) == rgn.getUpper(i)+ug[i]);
  }

  for (int i=elcField.getLowerExt(0); i<elcField.getUpperExt(0); ++i)
    for (int j=elcField.getLowerExt(1); j<elcField.getUpperExt(1); ++j)
      for (int k=elcField.getLower(2); k<elcField.getUpper(2); ++k)
        LC_ASSERT("Testing default values", elcField(i,j,k) == 10.0);

  Lucee::FieldPtr<double> ptr = elcField.createPtr();
  for (int i=elcField.getLowerExt(0); i<elcField.getUpperExt(0); ++i)
    for (int j=elcField.getLowerExt(1); j<elcField.getUpperExt(1); ++j)
    {
      elcField.setPtr(ptr, i, j);
      for (unsigned n=0; n<ptr.getNumComponents(); ++n)
        ptr[n] = (i+0.5*j)*n;
    }

  for (int i=elcField.getLowerExt(0); i<elcField.getUpperExt(0); ++i)
    for (int j=elcField.getLowerExt(1); j<elcField.getUpperExt(1); ++j)
      for (int k=elcField.getLower(2); k<elcField.getUpper(2); ++k)
        LC_ASSERT("Testing set values", elcField(i,j,k) == (i+0.5*j)*k);

// make a copy of elcField
  Lucee::Field<2, double> elcField2(elcField);
  LC_ASSERT("Testing if number of components is correct", 
    elcField2.getNumComponents() == 3);

  localRgn = elcField2.getRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      localRgn.getLower(i) == rgn.getLower(i));
    LC_ASSERT("Testing if index region is correct",
      localRgn.getUpper(i) == rgn.getUpper(i));
  }

  extRgn = elcField2.getExtRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      extRgn.getLower(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      extRgn.getUpper(i) == rgn.getUpper(i)+ug[i]);
  }

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      elcField2.getLowerExt(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      elcField2.getUpperExt(i) == rgn.getUpper(i)+ug[i]);
  }

  for (int i=elcField2.getLowerExt(0); i<elcField2.getUpperExt(0); ++i)
    for (int j=elcField2.getLowerExt(1); j<elcField2.getUpperExt(1); ++j)
      for (int k=elcField2.getLower(2); k<elcField2.getUpper(2); ++k)
        LC_ASSERT("Testing set values", elcField2(i,j,k) == (i+0.5*j)*k);

  for (int i=elcField2.getLowerExt(0); i<elcField2.getUpperExt(0); ++i)
    for (int j=elcField2.getLowerExt(1); j<elcField2.getUpperExt(1); ++j)
      for (int k=elcField2.getLower(2); k<elcField2.getUpper(2); ++k)
        elcField2(i,j,k) = 5.0*(i+0.5*j)*k;

  for (int i=elcField.getLowerExt(0); i<elcField.getUpperExt(0); ++i)
    for (int j=elcField.getLowerExt(1); j<elcField.getUpperExt(1); ++j)
      for (int k=elcField.getLower(2); k<elcField.getUpper(2); ++k)
        LC_ASSERT("Testing set values", elcField(i,j,k) == 5.0*(i+0.5*j)*k);
}

void
test_12()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  int lg[2] = {1, 3};
  int ug[2] = {2, 5};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> elcField(rgn, 3, lg, ug, 10.0);

  Lucee::FieldPtr<double> ptr = elcField.createPtr();
  for (int i=elcField.getLowerExt(0); i<elcField.getUpperExt(0); ++i)
    for (int j=elcField.getLowerExt(1); j<elcField.getUpperExt(1); ++j)
    {
      elcField.setPtr(ptr, i, j);
      for (unsigned n=0; n<ptr.getNumComponents(); ++n)
        ptr[n] = (i+0.5*j)*(n+1);
    }

  Lucee::Field<2, double> Ex = elcField.getSubCompView(0, 1);
  Lucee::Region<2, int> localRgn = Ex.getRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      localRgn.getLower(i) == rgn.getLower(i));
    LC_ASSERT("Testing if index region is correct",
      localRgn.getUpper(i) == rgn.getUpper(i));
  }

  Lucee::Region<2, int> extRgn = Ex.getExtRegion();
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      extRgn.getLower(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      extRgn.getUpper(i) == rgn.getUpper(i)+ug[i]);
  }

  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing if index region is correct",
      Ex.getLowerExt(i) == rgn.getLower(i)-lg[i]);
    LC_ASSERT("Testing if index region is correct",
      Ex.getUpperExt(i) == rgn.getUpper(i)+ug[i]);
  }

  LC_ASSERT("Testing if number of components is correct", Ex.getNumComponents() == 1);
  Lucee::FieldPtr<double> exPtr = Ex.createPtr();
  for (int i=Ex.getLowerExt(0); i<Ex.getUpperExt(0); ++i)
    for (int j=Ex.getLowerExt(1); j<Ex.getUpperExt(1); ++j)
    {
      Ex.setPtr(exPtr, i, j);
      for (unsigned n=0; n<exPtr.getNumComponents(); ++n)
        LC_ASSERT("Testing Ex values", exPtr[n] == (i+0.5*j));
    }
  for (int i=Ex.getLowerExt(0); i<Ex.getUpperExt(0); ++i)
    for (int j=Ex.getLowerExt(1); j<Ex.getUpperExt(1); ++j)
      LC_ASSERT("Testing Ex values", Ex(i,j,0) == (i+0.5*j));

  Lucee::Field<2, double> Eyz = elcField.getSubCompView(1, 3);
  LC_ASSERT("Checking number of components", Eyz.getNumComponents() == 2);
  Lucee::FieldPtr<double> eyzPtr = Eyz.createPtr();
  for (int i=Eyz.getLowerExt(0); i<Eyz.getUpperExt(0); ++i)
    for (int j=Eyz.getLowerExt(1); j<Eyz.getUpperExt(1); ++j)
    {
      Eyz.setPtr(eyzPtr, i, j);
      for (unsigned n=0; n<eyzPtr.getNumComponents(); ++n)
        LC_ASSERT("Testing Eyz values", eyzPtr[n] == (i+0.5*j)*(n+1+1));
    }
  LC_ASSERT("Testing component bounds of Eyz", Eyz.getLower(2) == 0);
  LC_ASSERT("Testing component bounds of Eyz", Eyz.getUpper(2) == 2);
  LC_ASSERT("Testing component bounds of Eyz", Eyz.getLowerExt(2) == 0);
  LC_ASSERT("Testing component bounds of Eyz", Eyz.getUpperExt(2) == 2);

  for (int i=Eyz.getLowerExt(0); i<Eyz.getUpperExt(0); ++i)
    for (int j=Eyz.getLowerExt(1); j<Eyz.getUpperExt(1); ++j)
      for (int k=Eyz.getLowerExt(2); k<Eyz.getUpperExt(2); ++k)
        LC_ASSERT("Testing Eyz values", Eyz(i,j,k) == (i+0.5*j)*(k+1+1));
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
  test_6();
  test_7();
  test_8();
  test_9();
  test_10();
  test_11();
  test_12();
  LC_END_TESTS;
}
