/**
 * @file	LcSOLPositivityDragCellUpdater.h
 *
 * @brief	Updater used to adjust energy at each cell by use of a drag term
 */

#ifndef LC_SOL_POSITIVITY_DRAG_CELL_UPDATER_H
#define LC_SOL_POSITIVITY_DRAG_CELL_UPDATER_H

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
  class SOLPositivityDragCellUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLPositivityDragCellUpdater();

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
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/** Matrix used to compute total number in a cell */
      Eigen::MatrixXd densityMatrix;
/**
 * When multiplied by the solution in a 3d cell, returns the integrated cell quantity
 * To get average, one should divide by the volume of the cell
 */
      Eigen::VectorXd mom0Vector;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

/**
 * Determines if two nodes have the same configuration space coordinates
 */
      bool sameConfigCoords(int srcIndex, int tarIndex, double dxMin, const Eigen::MatrixXd& nodeList);
  };
}

#endif // LC_SOL_POSITIVITY_DRAG_CELL_UPDATER_H
