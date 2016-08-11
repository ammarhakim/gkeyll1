/**
 * @file	LcDistFuncMomentCalc2D.h
 *
 * @brief	Updater to compute 2d moments of a 4d distribution function.
 * Currently only works for 0th moment
 */

#ifndef LC_DIST_FUNC_MOMENT_CALC_2D_H
#define LC_DIST_FUNC_MOMENT_CALC_2D_H

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
  class DistFuncMomentCalc2D : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC4 = Lucee::CDIM<4>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC2 = Lucee::CDIM<2>::N;

/** Create new modal DG solver in 4D */
      DistFuncMomentCalc2D();
/**  Destructor */
      virtual ~DistFuncMomentCalc2D();

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
/** Pointer to 4D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<4> *nodalBasis4d;
/** Pointer to 2D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Moment to compute */
      unsigned calcMom;
/** Space to store moment data for on a processor */
      Lucee::Field<2, double> *momentLocal;
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

#endif // LC_DIST_FUNC_MOMENT_CALC_2D_H
