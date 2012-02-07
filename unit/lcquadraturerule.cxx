/**
 * @file	lcquadraturerule.cxx
 *
 * @brief	Unit tests for Lucee::QuadratureRule class
 */

// lucee includes
#include <LcGaussianQuadRule.h>
#include <LcQuadratureRule.h>
#include <LcTest.h>

// std include
#include <cmath>

void
test_1()
{
  Lucee::GaussianQuadRule gauss2(2);

  LC_ASSERT("Testing Gaussian rule", gauss2.getNumNodes() == 2);

  Lucee::Vector<double> w(2), mu(2);
  gauss2.getOrdinatesAndWeights(mu, w);

  LC_ASSERT("Testing if weights sum to 2", epsCmp(w[0]+w[1], 2.0));
  LC_ASSERT("Testing ordinate", epsCmp(mu[0], -1/std::sqrt(3)));
  LC_ASSERT("Testing ordinate", epsCmp(mu[1], 1/std::sqrt(3)));

  Lucee::GaussianQuadRule gauss3(3);

  LC_ASSERT("Testing Gaussian rule", gauss3.getNumNodes() == 3);

  Lucee::Vector<double> w3(3), mu3(3);
  gauss3.getOrdinatesAndWeights(mu3, w3);

  double EPS = 1e-14;

  LC_ASSERT("Testing if weights sum to 2", diffCmp(w3[0]+w3[1]+w3[2], 2.0, EPS));

  LC_ASSERT("Testing ordinate", diffCmp(w3[0], 5.0/9.0, EPS));
  LC_ASSERT("Testing ordinate", diffCmp(w3[1], 8.0/9.0, EPS));
  LC_ASSERT("Testing ordinate", diffCmp(w3[2], 5.0/9.0, EPS));

  LC_ASSERT("Testing ordinate", diffCmp(mu3[0], -std::sqrt(3./5.), EPS));
  LC_ASSERT("Testing ordinate", diffCmp(mu3[1], 0.0, EPS));
  LC_ASSERT("Testing ordinate", diffCmp(mu3[2], std::sqrt(3./5.), EPS));
}

int
main(int argc, char **argv) 
{
  LC_BEGIN_TESTS("lcquadraturerule");
  test_1();
  LC_END_TESTS;
}
