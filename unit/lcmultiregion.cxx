/**
 * @file	lcmultiregion.cxx
 *
 * @brief	Unit tests for Lucee::MultiRegion class
 */

// lucee includes
#include <LcMultiRegion.h>
#include <LcTest.h>

void
test_1()
{
// create connections for 3 regions connected together as follows
//
//
// +------------+-------------+
// |            |             |
// |     0      |     1       |
// |            |             |
// +------------+-------------+
//              |             |
//              |     2       |
//              |             |
//              +-------------+

  int shape[2];
  Lucee::MultiRegion<2, int> multiRgn;
// add regions
  shape[0] = 20; shape[1] = 10;
  multiRgn.addRegion(0, Lucee::Region<2, int>(shape));

  shape[0] = 30; shape[1] = 10;
  multiRgn.addRegion(1, Lucee::Region<2, int>(shape));

  shape[0] = 30; shape[1] = 15;
  multiRgn.addRegion(2, Lucee::Region<2, int>(shape));

  LC_ASSERT("Checking number of regions", multiRgn.getNumRegions() == 3);

// loop over regions, ensuring that each is unconnected at this point
  Lucee::MultiRegion<2, int>::iterator itr(multiRgn);
  for ( ; ! itr.atEnd(); itr.next() )
  {
  }
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcmultiregion");
  test_1();
  LC_END_TESTS;
}
