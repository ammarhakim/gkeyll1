/**
 * @file	lceigen.cxx
 *
 * @brief	Unit tests for using EIGEN in Lucee
 */

// lucee includes
#include <LcTest.h>

// eigen inlcudes
#include <Eigen/Core>

void
test_1()
{
  Eigen::MatrixXd m(4,4);
  for (unsigned i=0; i<4; ++i)
    for (unsigned j=0; j<4; ++j)
      m(i,j) = 1.0;

  LC_ASSERT("Testing rows", m.rows() == 4);
  LC_ASSERT("Testing cols", m.cols() == 4);

  m.resize(3,5);
  LC_ASSERT("Testing rows", m.rows() == 3);
  LC_ASSERT("Testing cols", m.cols() == 5);
}

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lceigen");
  test_1();
  LC_END_TESTS;
}
