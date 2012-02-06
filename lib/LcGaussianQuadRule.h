/**
 * @file	LcGaussianQuadRule.h
 *
 * @brief	Guassian quadrature.
 */

#ifndef LC_GAUSSIAN_QUAD_RULE_H
#define LC_GAUSSIAN_QUAD_RULE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcQuadratureRule.h>

namespace Lucee
{
/**
 * Gaussian quadrature on interval [-1,1].
 */
  class GaussianQuadRule : public Lucee::QuadratureRule
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create new Gaussian quadrature object.
 */
      GaussianQuadRule();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Returns weights and ordinates of quadrature rule. The supplied
 * vectors must be pre-allocated.
 *
 * @param ord On output, ordinates for quadrature.
 * @param wth On output, weights for quadrature.
 */
      virtual void getOrdinatesAndWeights(Lucee::Vector<double>& ord,
        Lucee::Vector<double>& wth);

    private:
/** Weights */
      Lucee::Vector<double> weights;
/** Ordinates */
      Lucee::Vector<double> ordinates;
  };
}

#endif // LC_QUADRATURE_RULE_H
