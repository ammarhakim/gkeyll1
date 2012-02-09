/**
 * @file	LcModalDg1DUpdater.h
 *
 * @brief	Updater to solver 1D hyperbolic equations using modal DG scheme
 */

#ifndef LC_MODAL_DG_1D_UPDATER_H
#define LC_MODAL_DG_1D_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve 1D hyperbolic equations using modal DG scheme.
 */
  class ModalDg1DUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      ModalDg1DUpdater();

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
/** Equation to solve */
      Lucee::HyperEquation *equation;
/** CFL number to use */
      double cfl;
/** Maximum CFL number */
      double cflm;
/** Number of basis functions to use */
      double numBasis;
/** Values of Legendre polynomials at the ordinates */
      Matrix<double> Pmk; // m <- Polynomial order k <- ordinate index
/** Values of Legendre polynomials derivate at the ordinates */
      Matrix<double> DPmk; // m <- Polynomial order k <- ordinate index
/** Normalization coefficients */
      Vector<double> normCoeff;
/** Weights for quadrature */
      Lucee::Vector<double> w;
/** Ordinates for quadrature */
      Lucee::Vector<double> mu;
  };
}

#endif // LC_MODAL_DG_1D_UPDATER_H
