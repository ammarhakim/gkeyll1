/**
 * @file	lclincombiner.cxx
 *
 * @brief	Unit tests for LinCombiner updater
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcField.h>
#include <LcLinCombiner.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  int lg[2] = {1, 3};
  int ug[2] = {2, 5};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> a(rgn, 3, lg, ug, 10);
  Lucee::Field<2, double> b(rgn, 3, lg, ug, 10);
  Lucee::Field<2, double> c(rgn, 3, lg, ug, 10);

  Lucee::FieldPtr<double> aPtr = a.createPtr();
  Lucee::FieldPtr<double> bPtr = b.createPtr();
  for (int i=a.getLowerExt(0); i<a.getUpperExt(0); ++i)
    for (int j=a.getLowerExt(1); j<a.getUpperExt(1); ++j)
    {
      a.setPtr(aPtr, i, j);
      for (unsigned n=0; n<3; ++n)
        aPtr[n] = (i+10*j)*((double) n+3);
      
      b.setPtr(bPtr, i, j);
      for (unsigned n=0; n<3; ++n)
        bPtr[n] = 3*(i+10*j)*((double) n+3);
    }

  Lucee::LinCombiner<2> linCombiner;
// call initialization sequence
  linCombiner.declareTypes();
  linCombiner.initialize();

  std::vector<double> coeff(2);
  coeff[0] = 0.5; coeff[1] = 1.5;
  linCombiner.setCoeff(coeff);

  std::vector<const Lucee::DataStructIfc*> inpVars(2);
  std::vector<Lucee::DataStructIfc*> outVars(1);

// set input, output data-structures
  inpVars[0] = &a; inpVars[1] = &b;
  outVars[0] = &c;
  linCombiner.setInpVars(inpVars);
  linCombiner.setOutVars(outVars);

  linCombiner.setCurrTime(0.0);
// run updater
  linCombiner.update(1.0);

  Lucee::FieldPtr<double> cPtr = c.createPtr();
  for (int i=a.getLowerExt(0); i<a.getUpperExt(0); ++i)
    for (int j=a.getLowerExt(1); j<a.getUpperExt(1); ++j)
    {
      c.setPtr(cPtr, i, j);
      for (unsigned n=0; n<3; ++n)
      {
        double v1 = coeff[0]*(i+10*j)*((double) n+3);
        double v2 = coeff[1]*3*(i+10*j)*((double) n+3);
        LC_ASSERT("Testing linCombiner", cPtr[n] == v1+v2);
      }
    }
}

int
main(void)
{
  LC_BEGIN_TESTS("lclincombiner");
  test_1();
  LC_END_TESTS;
}
