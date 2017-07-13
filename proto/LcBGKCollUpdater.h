/**
 * @file	LcBGKCollUpdater.h
 *
 * @brief	Updater to calculate RHS BGK source term. Cross temperature and bulk velocity based on the Greene (1973, Phys. Ffuids)
 */

#ifndef LC_BGK_COLL_UPDATER_H
#define LC_BGK_COLL_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// math include
#define _USE_MATH_DEFINES
#include <math.h>

namespace Lucee
{
/**
 * Maxwellian distribution.
 */
  template <unsigned CDIM, unsigned VDIM>
  class BGKCollUpdater : public Lucee::UpdaterIfc
  {
    // Number of components for coordinate arrays etc.
    static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
    static const unsigned CNC = Lucee::CDIM<CDIM>::N;
    public:
      /** Class id: this is used by registration system */
      static const char *id;

      /** Constructor */
      BGKCollUpdater();
      /** Destructor */
      virtual ~BGKCollUpdater();
     
/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

/**
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      Lucee::UpdaterStatus update(double t);

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      void declareTypes();

      double evaluateMaxwell(double n, 
			     double v[VDIM], 
			     double invVt,
			     double invVt2);

    private:
     
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
      /** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;

      bool sameConfigCoords(unsigned n, unsigned cn, double dxMin,
        const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC);

      /** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> phaseConfMap;

      // loaded physical constants
      double mass, elemCharge, permitivity, floor;

      // Maxwell distribution normalization factor 1/sqrt(2*pi)
      double maxwellNorm;
      // collision frequency constant factor
      double collFreqConst;
      // plasma parameter constant factor
      double plasmaParamConst;
  };
}

#endif // LC_MAXWELL_DIST_INIT_H
