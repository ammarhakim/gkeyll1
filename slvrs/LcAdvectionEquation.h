/**
 * @file	LcAdvectionEquation.h
 *
 * @brief	Advection equation
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
#ifndef LC_ADVECTION_EQUATION_H
#define LC_ADVECTION_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>

namespace Lucee
{
/**
 * Represents an Advection equation.
 */
  class AdvectionEquation : public Lucee::HyperEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      AdvectionEquation();

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
      void flux(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f);

/**
 * Compute the wave speeds in the system.
 *
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, this constains the speeds.
 */
      void speeds(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& s);

/**
 * Decompose jump into waves and wave-speeds using right and left states.
 *
 * @param jump Jump to decompose.
 * @param ql Left state conserved variables.
 * @param qr Right state conserved variables.
 * @param waves On output, waves. This matrix has shape (meqn X mwave).
 * @param s On output, wave speeds.
 */
      virtual void waves(const Lucee::ConstFieldPtr<double>& jump,
        const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to primitive variables.
 * @param v On output, primitive variables.
 */
      void primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v);

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to conserved variables.
 * @param q On output, conserved variables.
 */
      void conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q);

    private:
/** advection speeds */
      double u[3];
  };
}

#endif //  LC_ADVECTION_EQUATION_H
