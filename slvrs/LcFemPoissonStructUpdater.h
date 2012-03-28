/**
 * @file	LcFemPoissonStructUpdater.h
 *
 * @brief	Updater to solve Poisson equations with FEM scheme on stuctured grids.
 */

#ifndef LC_FEM_POISSON_STRUCT_UPDATER_H
#define LC_FEM_POISSON_STRUCT_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// petsc includes
#include <petsc.h>

// std includes
#include <map>
#include <vector>

namespace Lucee
{
/**
 * Updater to solve Poisson equations using the finite element
 * method. This updater works on structured grids.
 */
  template <unsigned NDIM>
  class FemPoissonStructUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      FemPoissonStructUpdater();

/**
 * Destroy updater.
 */
      ~FemPoissonStructUpdater();

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
      Lucee::NodalFiniteElementIfc<NDIM> *nodalBasis;
/** Flag to indicate if nodes in source are shared or not */
      bool srcNodesShared;
/** Petsc matrices to store linear operator */
      Mat stiffMat;
/** Petsc vectors for source and initial guess */
      Vec globalSrc, initGuess;
/** Krylov subspace method context */
      KSP ksp;

/** Structure to store BC data. */
      struct FemPoissonBcData
      {
/** Create object to indicate Bc was not set */
          FemPoissonBcData()
            : isSet(false)
          {}

/** Flag to indicate if Bc was set */
          bool isSet;
/** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
          unsigned type;
/** Value to apply */
          double value;
      };

/** Boundary conditions on left/right edges in each direction */
      FemPoissonBcData bc[NDIM][2];

/** Map of rows to Dirichlet BC values */
      std::map<int, double> rowBcValues;

/**
 * Function to parse out BC.
 *
 * @param lt Lua table to get BC data from.
 * @return Boundary condition data.
 */
      FemPoissonBcData getBcData(const Lucee::LuaTable& lt) const;
  };
}

#endif // LC_FEM_POISSON_STRUCT_UPDATER_H
