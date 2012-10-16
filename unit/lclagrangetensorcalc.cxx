/**
 * @file	lclagrangetensorbasiscalc.cxx
 *
 * @brief	Unit tests for Lagrange tensor tests.
 */

// lucee includes
#include <LcLagrangeTensorBasisCalc.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::LagrangeTensorBasisCalc<1> basis;
  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == 1);

  unsigned numNodes[1] = {4}; // 4-node element
  basis.calc(Lucee::GAUSSIAN, numNodes);

  LC_ASSERT("Checking number of nodes", basis.getNumNodes() == 4);
}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lclagrangetensorcalc");
  test_1();
  LC_END_TESTS;
}
