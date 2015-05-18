/**
 * @file	LcDistFuncMomentCalc1DFrom3D.h
 *
 * @brief	Updater to compute 1d moments of a 3d distribution function
 */

#ifndef LC_DIST_FUNC_MOMENT_CALC_1D_FROM_3D_H
#define LC_DIST_FUNC_MOMENT_CALC_1D_FROM_3D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcCDIM.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to compute moments of distribution function f(x,v)
 */
  class DistFuncMomentCalc1DFrom3D : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC3 = Lucee::CDIM<3>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC1 = Lucee::CDIM<1>::N;

/** Create new modal DG solver in 1D */
      DistFuncMomentCalc1DFrom3D();

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
/** Pointer to 3D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 1D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Moment to compute */
      unsigned calcMom;
/** Direction to calcuate moment in */
      unsigned momDir;
/** Zeroth moment matrix */
      Eigen::MatrixXd mom0Matrix;
/** First moment matrix */
      Eigen::MatrixXd mom1Matrix;
/** Second moment matrix */
      Eigen::MatrixXd mom2Matrix;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_DIST_FUNC_MOMENT_CALC_1D_FROM_3D_H
