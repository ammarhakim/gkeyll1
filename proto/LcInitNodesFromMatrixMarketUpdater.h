/**
 * @file	LcInitNodesFromMatrixMarketUpdater.h
 *
 * @brief	Given an input file (.mtx), reads data into a field
 */

#ifndef LC_INIT_NODES_FROM_MATRIX_MARKET_UPDATER_H
#define LC_INIT_NODES_FROM_MATRIX_MARKET_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>
#include <Eigen/SparseCore>

namespace Lucee
{
/**
 * Updater to set a nodal FE field from supplied Lua function
 */
  template <unsigned NDIM>
  class InitNodesFromMatrixMarketUpdater : public Lucee::UpdaterIfc
  {
// Number of components for coordinate arrays etc.
      static const unsigned NC = Lucee::CDIM<NDIM>::N;

    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      InitNodesFromMatrixMarketUpdater();

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
/** Stores total number of degrees of freedom in system */
      int totalNodes;
/** Data structure to which data will be read into from filename */
      Eigen::SparseMatrix<double> initialValueMatrix;
/** Name of input file to use (must be in same directory) */
      std::string filename;
  };
}

#endif // LC_INIT_NODES_FROM_MATRIX_MARKET_UPDATER_H
