/**
 * @file	LcAdvectionEquation.h
 *
 * @brief	Advection equation for testing 4D basis functions.
 * Aligned rectangular coordinate system has been assumed.
 */
#ifndef LC_ADVECTION_EQUATION_TEMPLATED_H
#define LC_ADVECTION_EQUATION_TEMPLATED_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>
#include <LcAlignedRectCoordSys.h>

namespace Lucee
{
/**
 * Represents an Advection equation.
 */
  template <unsigned NDIM>
  class AdvectionEquationTemplated : public Lucee::HyperEquation
  {
/** Numerical flux to use */
      enum NumFluxType { CENTRAL_FLUX, UPWIND_FLUX };

    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      AdvectionEquationTemplated();

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
 * @param auxVars Auxillary variables needed to compute fluxes.
 * @param f On output, this contains the flux.
 */
      void flux(const Lucee::RectCoordSys& c, const double* q, 
        const std::vector<const double*>& auxVars, double* f);

/**
 * Compute the minimum and maximum wave speeds in the system. s[0] is
 * the minimum wave speed and s[1] is the maximum wave speed.
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @param s On output, s[0] is the minimum speed and s[1] the maximum speed.
 */
      void speeds(const Lucee::RectCoordSys& c, const double* q, double s[2]);

/**
 * Compute primitive variables given conserved variables.
 *
 * @param q Conserved variables for which to compute primitive variables.
 * @param v On output, primitive variables.
 */
      void primitive(const double* q, double* v) const;

/**
 * Compute conserved variables given primitive variables.
 *
 * @param v Primitive variables for which to compute conserved variables.
 * @param q On output, conserved variables.
 */
      void conserved(const double* v, double* q) const;

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
 * @param auxVarsl Left auxillary variables needed to compute waves.
 * @param auxVarsr Right auxillary variables needed to compute waves.
 * @param waves On output, waves. This matrix has shape (meqn X mwave).
 * @param s On output, wave speeds.
 */
      virtual void waves(const Lucee::RectCoordSys& c,
        const Lucee::ConstFieldPtr<double>& jump,
        const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
        const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
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
 * @param auxVarsl Left auxillary variables needed to compute fluxes.
 * @param auxVarsr Right auxillary variables needed to compute fluxes.
 * @param f On output, this contains the numerical flux.
 * @return Maximum wave speed from left/right state.
 */
      double numericalFlux(const Lucee::RectCoordSys& c,
        const double* ql, const double* qr, 
        const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
        double* f);

/**
 * Project given vector on left-eigenvectors of flux-jacobian. The
 * conserved variables are already in the local coordinate system and
 * hence for isotropic equations the coordinate system 'c' can be
 * ignored.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param q Conserved variables at which flux Jacobian is to be computed.
 * @param vec Vector to project.
 * @param coeff On output, the projection of 'vec' on left-eigenvectors.
 */
      void projectOnLeftEigenvectors(const Lucee::RectCoordSys& c,
        const double* q, const double* vec, double* coeff);

/**
 * Reconstruct vector by weighted sum of right eigenvectors. The
 * conserved variables are already in the local coordinate system and
 * hence for isotropic equations the coordinate system 'c' can be
 * ignored.
 *
 * @param c Coordinate system in which to compute eigensystem.
 * @param q Conserved variables at which flux Jacobian is to be computed.
 * @param coeff Coefficients to multiply corresponding right-eigenvectors.
 * @param vec On output, the reconstructured vector.
 */
      void reconWithRightEigenvectors(const Lucee::RectCoordSys& c,
        const double* q, const double* coeff, double* vec);

/**
 * Check if conserved variables satisfies invariant domains of the
 * system. Return true if it does, false otherwise.
 *
 * @param q Conserved variables.
 * @return true if invariant domains are satisfied, false otherwise.
 */
      bool isInvariantDomain(const double* q) const;

    private:
/** advection speeds */
      double u[NDIM];
/** Type of numerical flux to use */
      NumFluxType fluxType;
  };
}

#endif //  LC_ADVECTION_EQUATION_TEMPLATED_H
