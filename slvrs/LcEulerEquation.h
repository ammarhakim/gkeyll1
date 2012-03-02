/**
 * @file	LcEulerEquation.h
 *
 * @brief	Euler equations for gas-dynamics.
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
      virtual void flux(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f);

/**
 * Compute the minimum and maximum wave speeds in the system. s[0] is
 * the minimum wave speed and s[1] is the maximum wave speed.
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, s[0] is the minimum speed and s[1] the maximum speed.
 */
      virtual void speeds(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& q, double s[2]);

/**
 * Compute the absolute maximum wave speed
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @return maxium absolute speed.
 */
      double maxAbsSpeed(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& q);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to compute primitive variables.
 * @param v On output, primitive variables.
 */
      virtual void primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v) const;

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to compute conserved variables.
 * @param q On output, conserved variables.
 */
      virtual void conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q) const;

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

/**
 * Compute numerical flux for this equation system. Numerical flux
 * depends on left and right states. This method should also return
 * the maximum wave speed computed from the states. The states are
 * already in the normal-tangent space and so for isotropic systems
 * the coordinate systems can be ignored.
 *
 * @param c Coordinate system in which to compute flux.
 * @param ql Left conserved variable state.
 * @param qr Right conserved variable state.
 * @param f On output, this contains the numerical flux.
 * @return Maximum wave speed from left/right state.
 */
      double numericalFlux(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        Lucee::FieldPtr<double>& f);

    protected:

    private:
/** Gas gamma */
      double gas_gamma;
/** Should we enforce positivity? */
      bool correct;
/** Minium pressure in case correct=true */
      double minPressure;
/** Minium density in case correct=true */
      double minDensity;

/**
 * Compute pressure from conserved variables.
 *
 * @param q conserved variables.
 * @return pressure
 */
      double pressure(const Lucee::ConstFieldPtr<double>& q) const;

/**
 * Get density with basement fix.
 *
 * @param rho Density to correct.
 * @param correct density
 */
      double getSafeRho(double rho) const;
  };
}

#endif //  LC_EULER_EQUATION_H
