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
/** Class id: this is used by registration system */
      static const char *id;

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
 * Return number of equations in system.
 *
 * @return number of equations.
 */
      unsigned getNumEqns() const { return meqn; }

/**
 * Return number of waves in system.
 *
 * @return number of waves.
 */
      unsigned getNumWaves() const { return mwave; }

/**
 * Compute flux for this equation system.
 *
 * @param q Conserved variables for which to compute flux.
 * @param f On output, this contains the flux.
 */
      virtual void flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f);

    protected:

    private:
/** Number of equations */
      unsigned meqn;
/** Number of waves */
      unsigned mwave;
  };
}

#endif //  LC_HYPER_EQUATION_H
