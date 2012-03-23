/**
 * @file	LcNodalPoissonBracketUpdater.h
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

#ifndef LC_NODAL_POISSON_BRACKET_UPDATER_H
#define LC_NODAL_POISSON_BRACKET_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve 2D Poisson bracket operator PDE of the form
 *
 *    da
 *   ---- + {a,b} = 0
 *    dt
 *
 * where {a,b} is the Poisson bracket operator. This updater updates
 * the field 'a' using a first-order forward Euler step with the
 * Poisson bracket computed using a DG scheme.
 */
  class NodalPoissonBracketUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      NodalPoissonBracketUpdater();

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
      Lucee::NodalFiniteElementIfc<2> *nodalBasis;
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Differentiation matrix in X-direction */
      Lucee::Matrix<double> diffMatrix_x;
/** Differentiation matrix in X-direction */
      Lucee::Matrix<double> diffMatrix_y;
/** Stiffness in X-direction */
      Lucee::Matrix<double> stiffMatrix_x;
/** Stiffness in Y-direction */
      Lucee::Matrix<double> stiffMatrix_y;

/**
 * Compute gradient in x-direction.
 *
 * @param phiK potential at nodes
 * @param phiPrimeK On output, gradient in x-direction
 */
      void calcGradient_x(const std::vector<double>& phiK,
        std::vector<double>& phiPrimeK);

/**
 * Compute gradient in y-direction.
 *
 * @param phiK potential at nodes
 * @param phiPrimeK On output, gradient in x-direction
 */
      void calcGradient_y(const std::vector<double>& phiK,
        std::vector<double>& phiPrimeK);
  };
}

#endif // LC_NODAL_POISSON_BRACKET_UPDATER_H
