/**
 * @file        LcDivEquation.h
 *
 * @brief	Divergence of flux specifed as auxillary variable.
 */
#ifndef LC_DIV_AUX_EQUATION_H
#define LC_DIV_AUX_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>

namespace Lucee
{
/**
 * Represents an hyperbolic equation which is divergence of a flux
 * specified as an auxillary variable.
 */
  template <unsigned NDIM>
  class DivEquation : public Lucee::HyperEquation
  {
/** Type of flux to use */
      enum DiffFluxType { THREE_POINT, FIVE_POINT };

    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      DivEquation();

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
      virtual void flux(const Lucee::RectCoordSys& c, const double* q, 
        const std::vector<const double*>& auxVars, double* f);


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
      virtual double numericalFlux(const Lucee::RectCoordSys& c,
        const double* ql, const double* qr, 
        const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
        double* f);

    private:
/** Type of flux to use */
      DiffFluxType fluxType;
  };
}

#endif //  LC_DIV_AUX_EQUATION_H
