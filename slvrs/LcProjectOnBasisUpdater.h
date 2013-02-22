/**
 * @file	LcProjectOnBasisUpdater.h
 *
 * @brief	Project a function of a basis functions.
 */

#ifndef LC_PROJECT_ON_BASIS_UPDATER_H
#define LC_PROJECT_ON_BASIS_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Updater to project a supplied Lua function onto a set of basis
 * functions. The coefficients of the basis function are stored in
 * row-major order.
 */
  template <unsigned NDIM>
  class ProjectOnBasisUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;


/** Create new projection updater */
      ProjectOnBasisUpdater();

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
/** Number of basis functions to project on */
      unsigned numBasis;
/** Number of gaussian quadrature points to use */
      unsigned numGaussPoints;
/** Reference to function to project */
      int fnRef;
/** Values of Legendre polynomials at the ordinates */
      Matrix<double> Pmk; // m <- Polynomial order k <- ordinate index
/** Weights for quadrature */
      Lucee::Vector<double> w;
/** Ordinates for quadrature */
      Lucee::Vector<double> mu;

/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param loc Location at which to evaluate function.
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        const double loc[3], std::vector<double>& res);
  };
}

#endif // LC_PROJECT_ON_BASIS_UPDATER_H
