/**
 * @file	LcDistFuncReflectionBcUpdater.h
 *
 * @brief	Applies particle refection BCs to (electron) distribution function
 */

#ifndef LC_SOL_3D_ELECTROSTATIC_DIST_FUNC_REFLECTION_BC_UPDATER
#define LC_SOL_3D_ELECTROSTATIC_DIST_FUNC_REFLECTION_BC_UPDATER

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
  class SOL3DElectrostaticDistFuncReflectionBCUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOL3DElectrostaticDistFuncReflectionBCUpdater();

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
/** Pointer to 1D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Flag to indicate if cutoff velocities need to be computed */
      bool computeCutoffVelocities;
/** Flag to indicate if BCs should be applied to left edge */
      bool applyLeftEdge;
/** Flag to indicate if BCs should be applied to right edge */
      bool applyRightEdge;
/** Mapping for 180 degree rotations */
      std::vector<unsigned> rotMap;
/** Contains the right edge node numbers */
      std::vector<int> rightEdgeNodeNums;
/** Contains the left edge node numbers */
      std::vector<int> leftEdgeNodeNums;
/** Tolerance to which cutoff velocities should be found */
      double cutoffTolerance;
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussEdgeOrdinates;
      Eigen::MatrixXd gaussEdgeOrdinatesLeftEdge;
/** Weights for edge gaussian quadrature points */
      std::vector<double> gaussEdgeWeights;
      std::vector<double> gaussEdgeWeightsLeftEdge;
/** 
 * Interpolation matrix for bringing data that lives on the left or right edge
 * to quadrature points on same surface. It is constructed from the right edge
 * basis function evaluations but should work for the left edge as well.
 * The quadrature locations for these nodes are also in gaussSurfOrdinates and
 * gaussSurfWeights.
 */
      Eigen::MatrixXd edgeNodeInterpMatrix;
      Eigen::MatrixXd edgeNodeInterpMatrixLeftEdge;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_SOL_3D_ELECTROSTATIC_DIST_FUNC_REFLECTION_BC_UPDATER
