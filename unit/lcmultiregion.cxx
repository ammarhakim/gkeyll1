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
// test for connectivity class
  Lucee::MultiRegionConnectivity mrConn(0, 0, 0);

  LC_RAISES("Trying to reset to incorrect value", mrConn.reset(0, 0, 10), Lucee::Except);
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcmultiregion");
  test_1();
  LC_END_TESTS;
}
