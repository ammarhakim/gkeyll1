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
  int shape[2];
  Lucee::MultiRegionConnectivity lc[2], uc[2];

  Lucee::MultiRegion<2, int> multiRgn;
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

// add region 0
  shape[0] = 20; shape[1] = 10;
  lc[0] = Lucee::MultiRegionConnectivity(); // unconnected
  lc[1] = Lucee::MultiRegionConnectivity(); // unconnected

  uc[0] = Lucee::MultiRegionConnectivity(1, 0, Lucee::LOWER);
  uc[1] = Lucee::MultiRegionConnectivity(); // unconnected

  multiRgn.addRegion(0, Lucee::Region<2, int>(shape), lc, uc);

// add region 1
  shape[0] = 30; shape[1] = 10;
  lc[0] = Lucee::MultiRegionConnectivity(0, 0, Lucee::UPPER);
  lc[1] = Lucee::MultiRegionConnectivity(2, 1, Lucee::UPPER);

  uc[0] = Lucee::MultiRegionConnectivity(); // unconnected
  uc[1] = Lucee::MultiRegionConnectivity(); // unconnected

  multiRgn.addRegion(1, Lucee::Region<2, int>(shape), lc, uc);

// add region 2
  shape[0] = 30; shape[1] = 15;
  lc[0] = Lucee::MultiRegionConnectivity(); // unconnected
  lc[1] = Lucee::MultiRegionConnectivity(); // unconnected

  uc[0] = Lucee::MultiRegionConnectivity(1, 1, Lucee::LOWER);
  uc[1] = Lucee::MultiRegionConnectivity(); // unconnected

  multiRgn.addRegion(2, Lucee::Region<2, int>(shape), lc, uc);
    
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcmultiregion");
  test_1();
  LC_END_TESTS;
}
