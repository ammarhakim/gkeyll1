/**
 * @file	LcTenMomentEquation.h
 *
 * @brief	Ten-moment equations for gas-dynamics.
 */
#ifndef LC_TEN_MOMENT_EQUATION_H
#define LC_TEN_MOMENT_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>

namespace Lucee
{
/**
 * Represents an Ten-moment equation of gas dynamics.
 */
  class TenMomentEquation : public Lucee::HyperEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      TenMomentEquation();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Rotate data to local coordinate system.
 *
 * @param c Coordinate system to rotate data to.
 * @param inQ Input conserved variables.
 * @param outQ Rotated conserved variables. 
 */
      void rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ);

/**
 * Rotate data to global coordinate system.
 *
 * @param c Coordinate system to rotate data to.
 * @param inQ Input conserved variables.
 * @param outQ Rotated conserved variables. 
 */
      void rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ);

/**
 * Compute flux for this equation system.
 *
 * @param c Coordinate system in which to compute flux.
 * @param q Conserved variables for which to compute flux.
 * @param f On output, this contains the flux.
 */
      virtual void flux(const Lucee::RectCoordSys& c, const double* q, 
        const std::vector<const double*>& auxVars, double* f);

/**
 * Compute the minimum and maximum wave speeds in the system. s[0] is
 * the minimum wave speed and s[1] is the maximum wave speed.
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, s[0] is the minimum speed and s[1] the maximum speed.
 */
      virtual void speeds(const Lucee::RectCoordSys& c, const double* q, double s[2]);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to compute primitive variables.
 * @param v On output, primitive variables.
 */
      virtual void primitive(const double* q, double* v) const;

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to compute conserved variables.
 * @param q On output, conserved variables.
 */
      virtual void conserved(const double* v, double* q) const;

/**
 * Decompose jump into waves and wave-speeds using right and left
 * states. The states and jump are already in the local coordinate
 * system specified by 'c'. Hence, in most case (equation system is
 * isotropic) the coordinate system should be ignored.
 *
 * @param c Coordinate system in which to compute waves.
 * @param jump Jump to decompose.
 * @param ql Left state conserved variables.
 * @param qr Right state conserved variables.
 * @param waves On output, waves. This matrix has shape (meqn X mwave).
 * @param s On output, wave speeds.
 */
      virtual void waves(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& jump,
        const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s);

    private:
  };
}

#endif //  LC_EULER_EQUATION_H
