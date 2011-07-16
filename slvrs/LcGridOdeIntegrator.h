/**
 * @file   LcGridOdeIntegrator.h
 *
 * @brief   Base class for ODE integrator over complete grid.
 */

#ifndef LC_GRID_ODE_INTEGRATOR_H
#define LC_GRID_ODE_INTEGRATOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcField.h>
#include <LcStructuredGridBase.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Base class for ODE integrator for functions defined on grids.
 */
  template <unsigned NDIM>
  class GridOdeIntegrator : public Lucee::BasicObj
  {
    public:
/**
 * Create ODE integrator to integrate equations on specified grid.
 *
 * @param grid Grid to integrate ODEs on.
 * @param numInp Number of input fields to expect.
 */
      GridOdeIntegrator(const Lucee::StructuredGridBase<NDIM>& grid, unsigned numInp);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Integrate ODEs to time "t".
 *
 * @param t Time to integrate ODEs to.
 * @param sol On output, solution of ODEs.
 */
      virtual void integrate(double t, const Lucee::Field<NDIM, double>& sol) = 0;

/**
 * Get reference to 'loc' input field needed in ODE integration.
 *
 * @param loc Input field to fetch.
 * @return Reference to field.
 */
      const Lucee::Field<NDIM, double>& getIn(unsigned loc) const;

/**
 * Set field at location 'loc'.
 *
 * @param loc Location of field to set.
 * @param in Reference to field to set.
 */
      void setIn(unsigned loc, const Lucee::Field<NDIM, double>& in);

/**
 * Get reference to grid on which ODEs need to be solved.
 *
 * @return reference to grid.
 */
      const Lucee::StructuredGridBase<NDIM>& getGrid() const
      {
        return *gridPtr;
      }

    private:
/** Number of input fields */
      unsigned numInp;
/** Grid on which ODEs need to be solved */
      const Lucee::StructuredGridBase<NDIM> *gridPtr;
/** List of input fields */
      std::vector<const Lucee::Field<NDIM, double>* > inpFlds;
  };
}

#endif // LC_GRID_ODE_INTEGRATOR_H
