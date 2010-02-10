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
}

int
main(void)
{
  LC_BEGIN_TESTS("lcfield");
  test_1();
  LC_END_TESTS;
}
