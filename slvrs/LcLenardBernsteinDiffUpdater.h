/**
 * @file	LcLenardBernsteinDiffUpdater.h
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator.
 */

#ifndef LC_LENARD_BERNSTEIN_DIFF_UPDATER_H
#define LC_LENARD_BERNSTEIN_DIFF_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to evaluate the diffusion term in the Lenard-Bernstein
 * collision operator.
 */
  class LenardBernsteinDiffUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      LenardBernsteinDiffUpdater();

/**
 * Destroy updater.
 */
      ~LenardBernsteinDiffUpdater();

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
/** CFL number */
      double cfl;
/** Direction to apply diffusion operator */
      int diffDir;
/** Flag to indicate if to compute only increments */
      bool onlyIncrement;
/** Flag to indicate if Braginskii collisional time should be used */
      bool useBraginskii;
/** The following (in SI units) is used for Braginskii collisional time */
      double ionMass;
      double elementaryCharge;
      double epsilon0;
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** List of matrices for current cell */
      std::vector<Lucee::Matrix<double> > iMat;
/** List of matrices on each lower face */
      std::vector<Lucee::Matrix<double> > lowerMat;
/** List of matrices on each upper face */
      std::vector<Lucee::Matrix<double> > upperMat;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussSurfOrdinates;
/** Weights for surface gaussian quadrature points */
      std::vector<double> gaussSurfWeights;
/** 
 * Interpolation matrix for bringing data that lives on surface in diffDir
 * to quadrature points on same surface. It is constructed from the bottom surf
 * basis function evaluations but should work for upper surface as well.
 * The quadrature locations for these nodes are also in gaussSurfOrdinates and
 * gaussSurfWeights.
 */
      Eigen::MatrixXd surfNodeInterpMatrix;
/**
 * Compute matrix-vector multiply. Output vector must be
 * pre-allocated. Note that the computation performed is
 *
 * out = m*mat*vec + v*out
 *
 * @param m Factor for accumulation.
 * @param mat Matrix for multiplication.
 * @param vec Vector for multiplication.
 * @param v Factor for accumulation.
 * @param out On output, holds the product.
 */
      void matVec(double m, const Lucee::Matrix<double>& mat,
        const double* vec, double v, double* out);

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_LENARD_BERNSTEIN_DIFF_UPDATER_H
