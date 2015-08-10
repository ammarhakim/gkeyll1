/**
 * @file	LcNodalVlasovUpdater.h
 *
 * @brief	Updater to solve Vlasov equations with nodal DG scheme.
 */

#ifndef LC_NODAL_VLASOV_UPDATER_H
#define LC_NODAL_VLASOV_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve Vlasov equation using a nodal discontinous
 * Galerkin scheme.
 */
  template <unsigned CDIM, unsigned VDIM>
  class NodalVlasovUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      NodalVlasovUpdater();

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
/** Directions to update */
      std::vector<unsigned> updateDims;
/** Pointer to phase-space basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
/** Pointer to configuration-space basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Flag to indicate if to only compute increments */
      bool onlyIncrement;

/**
 * Matrix holder: this class is needed as the Matrix class does not
 * have a default constructor.
 */
      struct MatrixHolder
      {
/** Ctor */
          MatrixHolder() : m(1, 1) {}
/** Matrix data */
          Lucee::Matrix<double> m;
      };

/** Stiffness matrices */
      MatrixHolder stiffMatrix[CDIM+VDIM];
/** Lifting matrices for lower surface */
      MatrixHolder lowerLift[CDIM+VDIM];
/** Lifting matrices for upper surface */
      MatrixHolder upperLift[CDIM+VDIM];

/**
 * Structure to store node numbers on edges.
 */
      struct EdgeNodeNums
      {
/** Node numbers */
          std::vector<int> nums;
      };

/** Vector to store lower node numbers */
      EdgeNodeNums lowerNodeNums[CDIM+VDIM];
/** Vector to store upper node numbers */
      EdgeNodeNums upperNodeNums[CDIM+VDIM];

/**
 * Compute matrix-vector multiply. Output vector must be
 * pre-allocated. Note that the computation performed is
 *
 * out = m*mat*vec + v*out
 *
 * @param m Factor for accumulation.
 * @param mat Matrix for multiplication.
 * @param meqn Number of equations.
 * @param vec Vector for multiplication.
 * @param v Factor for accumulation.
 * @param out On output, holds the product.
 */
      void matVec(double m, const Lucee::Matrix<double>& mat,
        unsigned meqn, const double* vec, double v, double* out);
  };
}

#endif // LC_NODAL_VLASOV_UPDATER_H
