/**
 * @file	lcrectcartgrid.cxx
 *
 * @brief	Unit tests for Lucee::FixedVector class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcTest.h>
#include <LcRectCartGrid.h>
#include <LcRowMajorSequencer.h>

// std includes
#include <cmath>

void
test_1()
{
  int cells[3] = {1, 2, 4};
  Lucee::Region<3, int> localBox(cells);
  double lower[3] = {0.0, 0.0, 0.0};
  double upper[3] = {1.0, 2.0, 4.0};
  Lucee::Region<3, double> physBox(lower, upper);
  Lucee::RectCartGrid<3> grid3(localBox, localBox, physBox);

// cell volume
  int idx[3];
  Lucee::RowMajorSequencer<3> seq(localBox);
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    LC_ASSERT("Testing cell volume", grid3.getVolume() == 1.0);
  }

  LC_ASSERT("Testing number of cells", grid3.getNumCells(0) == 1);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(1) == 2);
  LC_ASSERT("Testing number of cells", grid3.getNumCells(2) == 4);

  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid3.setIndex(idx);
    for (unsigned i=0; i<3; ++i)
      LC_ASSERT("Testing surface area volume", grid3.getSurfArea(i) == 1.0);
  }
}

int
main(void)
{
  LC_BEGIN_TESTS("lcrectcartgrid");
  test_1();
  LC_END_TESTS;
}
