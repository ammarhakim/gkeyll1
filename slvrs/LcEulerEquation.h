/**
 * @file	LcEulerEquation.h
 *
 * @brief	Euler equations for gas-dynamics.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
#ifndef LC_EULER_EQUATION_H
#define LC_EULER_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>

namespace Lucee
{
/**
 * Represents an Euler equation of gas dynamics.
 */
  class EulerEquation : public Lucee::HyperEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      EulerEquation();

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
      virtual void flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f);

/**
 * Compute the wave speeds in the system.
 *
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, this constains the speeds.
 */
      virtual void speeds(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double> s);

    protected:

    private:
/** Gas gamma */
      double gas_gamma;
  };
}

#endif //  LC_EULER_EQUATION_H
