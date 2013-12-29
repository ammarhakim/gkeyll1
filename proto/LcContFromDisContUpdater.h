/**
 * @file	LcContFromDisContUpdater.h
 *
 * @brief	Updater to get a continous function from a discontinous one.
 */

#ifndef LC_CONT_FROM_DIS_CONT_UPDATER_H
#define LC_CONT_FROM_DIS_CONT_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// petsc includes
#include <petsc.h>

// std includes
#include <map>
#include <vector>

namespace Lucee
{
/**
 * Updater to project a discontinuous function on a continous set of
 * basis. This is a copy of the FemPoissonStructUpdater and hence is
 * templated over NDIM. However, the instantiation is just for
 * NDIM=1. The code is otherwise almost identical to the Poisson
 * solver..
 */
  template <unsigned NDIM>
  class ContFromDisContUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      ContFromDisContUpdater();

/**
 * Destroy updater.
 */
      ~ContFromDisContUpdater();

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

/** Map of rows to Dirichlet BC values */
      std::map<int, double> rowBcValues;

/** Flags to indicated periodic directions */
      bool periodicFlgs[NDIM];
/** Flag to indicate if all directions are periodic */
      bool allPeriodic;
/** Flag to indicate if update() method was called at least once */
      bool runOnce;
  };
}

#endif // LC_CONT_FROM_DIS_CONT_UPDATER_H
