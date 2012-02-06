/**
 * @file	LcQuadratureRule.h
 *
 * @brief	Base class for a quadrature rule.
 */

#ifndef LC_QUADRATURE_RULE_H
#define LC_QUADRATURE_RULE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Interface class for quadrature rule providing interface to get
 * weights and abscissa.
 */
  class QuadratureRule : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/** Destroy dataStruct */
      virtual ~QuadratureRule();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Return number of nodes in quadrature.
 *
 * @return Number of nodes in quadrature.
 */
      virtual void getNumNodes() const = 0;

/**
 * Returns weights and ordinates of quadrature rule. The supplied
 * vectors must be pre-allocated.
 *
 * @param ord On output, ordinates for quadrature.
 * @param wth On output, weights for quadrature.
 */
      virtual void getOrdinatesAndWeights(Lucee::Vector<double>& ord,
        Lucee::Vector<double>& wth) = 0;
  };
}

#endif // LC_QUADRATURE_RULE_H
