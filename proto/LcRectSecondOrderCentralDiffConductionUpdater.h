/**
 * @file	LcRectSecondOrderCentralDiffConductionUpdater.h
 *
 * @brief	Compute 2nd order central-differences on a rectangular grid.for anisotropic heat conductivity with B field
 */

#ifndef LC_RECT_SECOND_ORDER_CENTRAL_DIFF_CONDUCTION_UPDATER_H
#define LC_RECT_SECOND_ORDER_CENTRAL_DIFF_CONDUCTION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to compute second-order central difference of the temperature
 * defined on a rectangular grid decomposed into parallel and perpendicular
 * directions, scaled by a diffusion constant.
 */
  template <unsigned NDIM>
  class RectSecondOrderCentralDiffConductionUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new central-difference updater.
 */
      RectSecondOrderCentralDiffConductionUpdater();
      
/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

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

    private:
/**
 * Compute central heat flux using parallel and perpendicular conduction.
 *
 * @param inFld Field to compute CD of, in this case the pressure
 * @param rhoFld Fluid density
 * @param emFld Electromagnetic field
 * @param outFld Difference output 
 * @param qFld Heat flux tensor, used for calculating output
 */
      void computeConduction1D(const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld);


/**
 * Compute central heat flux using parallel and perpendicular conduction.
 *
 * @param inFld Field to compute CD of, in this case the pressure
 * @param rhoFld Fluid density
 * @param emFld Electromagnetic field
 * @param outFld Difference output 
 * @param qFld Heat flux tensor, used for calculating output
 *
 */
      void computeConduction2D(const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld);

/**
 * Compute central heat flux using parallel and perpendicular conduction
 *
 * @param inFld Field to compute CD of, in this case the pressure
 * @param rhoFld Fluid density
 * @param emFld Electromagnetic field
 * @param outFld Difference output 
 * @param qFld Heat flux tensor, used for calculating output
 */
      void computeConduction3D(const Lucee::Field<NDIM, double>& inFld, const Lucee::Field<NDIM, double>& rhoFld, const Lucee::Field<NDIM, double>& emField, Lucee::Field<NDIM, double>& outFld, Lucee::Field<NDIM, double>& qFld);

/**
 * Compute the symmetric partial P/rho  at a single point given temperature gradients and pressure. Scale this by the parallel and perp conductivities. 
 * @param p Pressure field
 * @param dT1 Temperature gradient divided by mass in direction 1
 * @param dT2 Temperature gradient divided by mass in direction 2
 * @param dT3 Temperature gradient divided by mass in direction 3
 * @param B Magnetic field vector
 * @param rho Density
 * @param q Output heat flux tensor (10 components)
 *
 */
      void heatFlux(const double* p, const double* dT1, const double* dT2, const double* dT3, const double* B, const double rho, double* q);


/** Cell spacing */
      double dx[NDIM];
/** scale factor for larmor radius in denominator of q_perp  */
      double rhoFactor;
/** species charge */
      double charge;
/** species mass */
      double mass;
/** conductivity/v_thermal */
      double alpha;

  };
}

#endif // LC_RECT_SECOND_ORDER_CENTRAL_DIFF_CONDUCTION_UPDATER_H
