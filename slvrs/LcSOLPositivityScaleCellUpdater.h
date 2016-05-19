/**
 * @file	LcSOLPositivityScaleCellUpdater.h
 *
 * @brief	Updater used to adjust energy at each cell by use of a drag term
 * Uses result of LcSOLPositivityDragCellUpdater to add the correct portion to match
 * the desired energy at each node.
 */

#ifndef LC_SOL_POSITIVITY_SCALE_CELL_UPDATER_H
#define LC_SOL_POSITIVITY_SCALE_CELL_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

// eigen includes
#include <Eigen/Core>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOLPositivityScaleCellUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLPositivityScaleCellUpdater();

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
/** Pointer to phase space basis functions to use */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Default true. If true, limits will be set on scale factors and positivity checks will be made */
      bool positivityChecks;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_SOL_POSITIVITY_SCALE_CELL_UPDATER_H
