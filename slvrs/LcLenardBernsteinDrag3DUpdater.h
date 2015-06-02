/**
 * @file	LcLenardBernsteinDrag3DUpdater.h
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 * df/dt = alpha*d[(v-u)f]/dv where u is a function of x (drift velocity)
 * Used for 1D2V SOL problem.
 */

#ifndef LC_LENARD_BERNSTEIN_DRAG_3D_UPDATER_H
#define LC_LENARD_BERNSTEIN_DRAG_3D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using a nodal discontinous
 * Galerkin scheme.
 */
  class LenardBernsteinDrag3DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      LenardBernsteinDrag3DUpdater();

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
/** Diffusion coefficient */
      double alpha;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis;
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Flag to indicate if to only compute increments */
      bool onlyIncrement;
/** Reference to function specifying alpha */
      int fnRef;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      std::vector<Eigen::MatrixXd> gaussSurfOrdinates;
/** Weights for surface gaussian quadrature points */
      std::vector<std::vector<double> > gaussSurfWeights;
/**
 * Matrix of volume gaussian quadrature locations.
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing volume integrals.
 */
      Eigen::MatrixXd gaussVolOrdinates;
/** Weights for volume gaussian quadrature points */
      std::vector<double> gaussVolWeights;
/** 
 * Derivative of basis functions in drag direction evaluated at
 * volume gaussian quadrature points. Rows correspond to different
 * quadrature points specified in gaussVolOrdinates and columns 
 * correspond to different basis functions
 */
      std::vector<Eigen::MatrixXd> basisDerivAtVolQuad;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to volume
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpVolMatrix;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to volume
 * gaussian quadrature points for a 1d field.
 */
      Eigen::MatrixXd interpVolMatrix1d;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to surface
 * gaussian quadrature points on lower face in dragDir.
 */
      std::vector<Eigen::MatrixXd> interpSurfMatrixLower;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to surface
 * gaussian quadrature points on upper face in dragDir.
 */
      std::vector<Eigen::MatrixXd> interpSurfMatrixUpper;
/** 
 * Mass matrix inverse times tranpose of interpSurfMatrixLower precomputed.
 */
      std::vector<Eigen::MatrixXd> surfIntegralMatrixLower;
/** 
 * Mass matrix inverse times tranpose of interpSurfMatrixUpper precomputed.
 */
      std::vector<Eigen::MatrixXd> surfIntegralMatrixUpper;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        std::vector<double>& res);
  };
}

#endif // LC_LENARD_BERNSTEIN_DRAG_3D_UPDATER_H
