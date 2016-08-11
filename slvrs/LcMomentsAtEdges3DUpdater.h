/**
 * @file	LcMomentsAtEdges3DUpdater.h
 *
 * @brief	Computes several parallel velocity moments of the distribution function at both edges.
 * Used for 3D SOL problem (Adiabatic Electrons).
 */

#ifndef LC_MOMENTS_AT_EDGES_3D_UPDATER_H
#define LC_MOMENTS_AT_EDGES_3D_UPDATER_H

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
#include <LcDynVector.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class MomentsAtEdges3DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      MomentsAtEdges3DUpdater();

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
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis;
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussEdgeOrdinatesLower;
      Eigen::MatrixXd gaussEdgeOrdinatesUpper;
/** Weights for edge gaussian quadrature points */
      std::vector<double> gaussEdgeWeightsLower;
      std::vector<double> gaussEdgeWeightsUpper;
/** Surface interpolation matrices */
      Eigen::MatrixXd interpEdgeMatrixLower;
      Eigen::MatrixXd interpEdgeMatrixUpper;
/**
 * Multiplies vector of nodal values to get derivatives at nodes
 */
      Eigen::MatrixXd derivativeMatrix;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_MOMENTS_AT_EDGES_3D_UPDATER_H
