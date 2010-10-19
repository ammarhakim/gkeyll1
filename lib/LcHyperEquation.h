/**
 * @file	LcHyperEquation.h
 *
 * @brief	Interface to hyperbolic equations.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
#ifndef LC_HYPER_EQUATION_H
#define LC_HYPER_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcConstFieldPtr.h>
#include <LcFieldPtr.h>

namespace Lucee
{
/**
 * Represents a hyperbolic equation.
 */
  class HyperEquation : public Lucee::BasicObj
  {
    public:
/**
 * Create a new hyperbolic equation system.
 *
 * @param Number of equations in system.
 * @param mwave Number of waves in system.
 */
      HyperEquation(unsigned meqn, unsigned mwave);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Compute flux for this equation system.
 *
 * @param q Conserved variables for which to compute flux.
 * @param f On output, this contains the flux.
 */
      virtual void flux(Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f);

    protected:

    private:
/** Number of equations */
      unsigned meqns;
/** Number of waves */
      unsigned mwave;
  };
}

#endif //  LC_HYPER_EQUATION_H
