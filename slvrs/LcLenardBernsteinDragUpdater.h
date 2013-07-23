/**
 * @file	LcLenardBernsteinDragUpdater.h
 *
 * @brief	Updater to compute the drag term in the L-B collision operator.
 * df/dt = alpha*d[(v-u)f]/dv where u is a function of x (drift velocity)
 * A key assumption I am making is that the basis functions in all the cells are the same.
 */

#ifndef LC_LENARD_BERNSTEIN_DRAG_UPDATER_H
#define LC_LENARD_BERNSTEIN_DRAG_UPDATER_H

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
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using a nodal discontinous
 * Galerkin scheme.
 */
  class LenardBernsteinDragUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      LenardBernsteinDragUpdater();

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
/** Directions to update */
      std::vector<unsigned> updateDims;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Direction to apply drag operator */
      int dragDir;
/** Flag to indicate if to only compute increments */
      bool onlyIncrement;
/** Flag to indicate if Braginskii collisional time should be used */
      bool useBraginskii;
/** The following (in SI units) is used for Braginskii collisional time */
      double ionMass;
      double elementaryCharge;
      double epsilon0;
/** Type of interface flux to use */
      int fluxType;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussSurfOrdinates;
/** Weights for surface gaussian quadrature points */
      std::vector<double> gaussSurfWeights;
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
      Eigen::MatrixXd basisDerivAtVolQuad;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to volume
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpVolMatrix;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to surface
 * gaussian quadrature points on lower face in dragDir.
 */
      Eigen::MatrixXd interpSurfMatrixLower;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to surface
 * gaussian quadrature points on upper face in dragDir.
 */
      Eigen::MatrixXd interpSurfMatrixUpper;
/** 
 * Mass matrix inverse times tranpose of interpSurfMatrixLower precomputed.
 */
      Eigen::MatrixXd surfIntegralMatrixLower;
/** 
 * Mass matrix inverse times tranpose of interpSurfMatrixUpper precomputed.
 */
      Eigen::MatrixXd surfIntegralMatrixUpper;
/** 
 * Interpolation matrix for bringing data that lives on surface in dragDir
 * to quadrature points on same surface. It is constructed from the bottom surf
 * basis function evaluations but should work for upper surface as well.
 * The quadrature locations for these nodes are also in gaussSurfOrdinates and
 * gaussSurfWeights.
 */
      Eigen::MatrixXd surfNodeInterpMatrix;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_LENARD_BERNSTEIN_DRAG_UPDATER_H
