/**
 * @file	LcGaussianQuadRule.cpp
 *
 * @brief	Guassian quadrature.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGaussianQuadRule.h>
#include <LcMathLib.h>

namespace Lucee
{
  const char *GaussianQuadRule::id = "Gaussian";

  GaussianQuadRule::GaussianQuadRule()
    : QuadratureRule(1), weights(1), ordinates(1)
  {
// set some reasonable defaults
    weights[0] = 2.0;
    ordinates[0] = 0.0;
  }

  GaussianQuadRule::GaussianQuadRule(unsigned n)
    : QuadratureRule(n), weights(n), ordinates(n)
  {
// compute weights and ordinates
    Lucee::gauleg(n, -1, 1, ordinates, weights);
  }

  void
  GaussianQuadRule::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    QuadratureRule::readInput(tbl);

// allocate space to compute weights and ordinates
    unsigned n = this->getNumNodes();
    weights = Lucee::Vector<double>(n);
    ordinates = Lucee::Vector<double>(n);

// compute weights and ordinates
    Lucee::gauleg(n, -1, 1, ordinates, weights);
  }

  void
  GaussianQuadRule::getOrdinatesAndWeights(Lucee::Vector<double>& ord, Lucee::Vector<double>& wth)
  {
    for (unsigned i=0; i<this->getNumNodes(); ++i)
    {
      ord[i] = ordinates[i];
      wth[i] = weights[i];
    }
  }
}
