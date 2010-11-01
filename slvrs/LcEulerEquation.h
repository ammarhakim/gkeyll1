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
      virtual void speeds(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& s);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to primitive variables.
 * @param v On output, primitive variables.
 */
      virtual void primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v);

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to conserved variables.
 * @param q On output, conserved variables.
 */
      virtual void conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q);

    protected:

    private:
/** Gas gamma */
      double gas_gamma;
/**
 * Compute pressure from conserved variables.
 *
 * @param q conserved variables.
 * @return pressure
 */
      double pressure(const Lucee::ConstFieldPtr<double>& q) const;
  };
}

#endif //  LC_EULER_EQUATION_H
