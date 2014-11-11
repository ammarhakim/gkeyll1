/**
 * @file	LcNodalCopy3DTo5DFieldUpdater.h
 *
 * @brief	Updater to copy data from a 3D field onto a 5D field using the same basis functions.
 * Currently only works assuming 3D field lies on the (0,1,2) surface
 * Currently hard-coded to assume Serendipity element.
 */

#ifndef LC_NODAL_COPY_3D_TO_5D_FIELD_UPDATER_H
#define LC_NODAL_COPY_3D_TO_5D_FIELD_UPDATER_H

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
 * Updater to copy a 2D nodal field to a 4D nodal field
 */
  class NodalCopy3DTo5DFieldUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

      NodalCopy3DTo5DFieldUpdater();

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
/** Pointer to 5d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to 3d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Matrix that maps the 3D nodes to 5D */
      Eigen::MatrixXd mappingMatrix;
/** Temporary flag to keep track of what polynomial order of element we are copying */
      int polyOrder;
  };
}

#endif // LC_NODAL_COPY_3D_TO_5D_FIELD_UPDATER_H
