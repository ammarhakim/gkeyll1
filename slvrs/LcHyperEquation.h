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
#include <LcMatrix.h>

namespace Lucee
{
/**
 * Represents a system of hyperbolic equations. This class provides a
 * rich interface to compute various quantities (fluxes, speeds,
 * waves, eigensystem) for use in various numercial schemes. Not all
 * methods are required for every scheme and the documentation for a
 * particular scheme should be consulted to determine which methods
 * are actually used.
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

/**
 * Compute the minimum and maximum wave speeds in the system. s[0] is
 * the minimum wave speed and s[1] is the maximum wave speed.
 *
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, s[0] is the minimum speed and s[1] the maximum speed.
 */
      virtual void speeds(const Lucee::ConstFieldPtr<double>& q, double s[2]);

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
 * Compute eigensystem for equations give a state. This method should
 * return all eigenvalues and right and left eigenvectors.
 *
 * @param q State at which to compute eigensystem.
 * @param ev On output, eigenvalues of system.
 * @param rev On output, right eigenvectors of system stored in columns.
 * @param lev On output, left eigenvectors of system stored in columns.
 */
      virtual void eigensystem(const Lucee::ConstFieldPtr<double>& q,
        Lucee::Vector<double>& ev, Lucee::Matrix<double>& rev, Lucee::Matrix<double>& lev);

/**
 * Check if conserved variables satisfies invariant domains of the
 * system. Return true if it does, false otherwise.
 *
 * @param q Conserved variables.
 * @return true if invariant domains are satisfied, false otherwise.
 */
      virtual bool isInvariantDomain(const Lucee::ConstFieldPtr<double>& q) const;

/**
 * Compute fluctuations using q-waves from waves and speeds.
 *
 * @param waves Waves. This matrix has shape (meqn X mwave).
 * @param s Wave speeds.
 * @param amdq On output, fluctuations in the negative direction.
 * @param apdq On output, fluctuations in the positive direction.
 */
      void qFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
        Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq);

/**
 * Compute fluctuations using f-waves from waves and speeds.
 *
 * @param waves Waves. This matrix has shape (meqn X mwave).
 * @param s Wave speeds.
 * @param amdq On output, fluctuations in the negative direction.
 * @param apdq On output, fluctuations in the positive direction.
 */
      void fFluctuations(const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
        Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq);

    protected:

    private:
/** Number of equations */
      unsigned meqn;
/** Number of waves */
      unsigned mwave;
  };
}

#endif //  LC_HYPER_EQUATION_H
