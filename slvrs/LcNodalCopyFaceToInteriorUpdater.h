/**
 * @file	LcNodalCopyFaceToInteriorUpdater.h
 *
 * @brief	Updater to copy data on the node face into the volume nodes.
 * Assumes common nodes are not shared
 */

#ifndef LC_NODAL_COPY_FACE_TO_INTERIOR_UPDATER_H
#define LC_NODAL_COPY_FACE_TO_INTERIOR_UPDATER_H

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
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to copy a 1D nodal field to a 2D nodal field
 */
  class NodalCopyFaceToInteriorUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

      NodalCopyFaceToInteriorUpdater();

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
/** Pointer to 2d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Pointer to 1d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Needed to distinguish continuous fields vs. discontinuous for now */
      bool shareCommonNodes;
/** Direction that 1-D nodal data "lives on" */
      unsigned dir;
/** Matrix that maps the 1-D nodes to 2-D */
      Eigen::MatrixXd mappingMatrix;
  };
}

#endif // LC_NODAL_COPY_FACE_TO_INTERIOR_UPDATER_H
