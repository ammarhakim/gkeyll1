/**
 * @file	LcNodalCopy1DTo3DFieldUpdater.h
 *
 * @brief	Updater to copy data from a 1D field onto a 3D field using the same basis functions.
 * Currently only works assuming 1D field lies along the 0 direction
 * Currently hard-coded to assume Serendipity element
 */

#ifndef LC_NODAL_COPY_1D_TO_3D_FIELD_UPDATER_H
#define LC_NODAL_COPY_1D_TO_3D_FIELD_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>

namespace Lucee
{
/**
 * Updater to copy a 1D nodal field to a 3D nodal field
 */
  class NodalCopy1DTo3DFieldUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

      NodalCopy1DTo3DFieldUpdater();

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

    private:
/** Pointer to 3d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 1d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Matrix that maps the 1D nodes to 3D */
      Eigen::MatrixXd mappingMatrix;
/** Temporary flag to keep track of what polynomial order of element we are copying */
      int polyOrder;
  };
}

#endif // LC_NODAL_COPY_1D_TO_3D_FIELD_UPDATER_H
