/**
 * @file	LcPhMaxwellEquation.h
 *
 * @brief	Perfectly-Hyperbolic Maxwell equations for electromagnetism.
 */
#ifndef LC_PH_MAXWELL_EQUATION_H
#define LC_PH_MAXWELL_EQUATION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcHyperEquation.h>

namespace Lucee
{
/**
 * Represents Perfectly-Hyperbolic Maxwell equations of
 * electromagnetism. This system explicity takes into account the
 * divergence relations.
 */
  class PhMaxwellEquation : public Lucee::HyperEquation
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new hyperbolic equation system.
 */
      PhMaxwellEquation();

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
 * Compute the absolute maximum wave speed
 *
 * @param c Coordinate system in which to compute speeds.
 * @param q Conserved variables for which to compute speeds.
 * @return maxium absolute speed.
 */
      double maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q);


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

/**
 * Check if conserved variables satisfies invariant domains of the
 * system. Return true if it does, false otherwise.
 *
 * @param q Conserved variables.
 * @return true if invariant domains are satisfied, false otherwise.
 */
      virtual bool isInvariantDomain(const double* q) const;

    protected:

    private:
/** Speed of light */
      double lightSpeed;
/** Speed of light squared */
      double lightSpeed2;
/** Propagation speed factor for electric field error potential. This
 * is dimensionless and the actual speed is chi_e*c. */
      double chi_e;
/** Propagation speed factor for magnetic field error potential. This
 * is dimensionless and the actual speed is chi_m*c. */
      double chi_m;
/** Maximum wave speed */
      double maxWaveSpeed;

/** Enum for flux types */
      enum NumFlux { NF_CENTRAL, NF_LAX, NF_UPWIND };

/** Flag to indicate type of numerical flux to use */
      NumFlux numFlux;
/** Flag to indicate use of intermediate wave (makes sense only if using Lax fluxes) */
      bool useIntermediateWave;
/** Flag to indicate if static EM field is present */
      bool hasStatic;
      
  };
}

#endif //  LC_PH_MAXWELL_EQUATION_H
