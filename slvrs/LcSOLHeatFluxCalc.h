/**
 * @file	LcSOLHeatFluxCalc.h
 *
 * @brief	Computes several parallel velocity moments of the distribution function at both edges.
 * Used for 5D SOL problem.
 */

#ifndef LC_SOL_HEAT_FLUX_CALC_H
#define LC_SOL_HEAT_FLUX_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcVector.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class SOLHeatFluxCalc : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      SOLHeatFluxCalc();

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
/** Pointer to phase space basis functions to use */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 2d (v,mu) space basis functions */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/** Keeps track of the offsets needed to get all nodes that share the same config. space location */
      std::vector<int> nodalStencil;
/** Flag indicating whether or not to integrate in ghost cells */
      bool integrateGhosts;
/**
 * Interpolation matrix used to compute integration in (v,mu) space
 */
      Eigen::MatrixXd interpMatrix2d;
/**
 * Quadrature weights used to compute integration in (v,mu) space
 */
      std::vector<double> gaussWeights2d;
/**
 *  Keeps track of nodes to write data to on lower surface in z
 */
      std::vector<int> lowerEdgeNodeNums;
/**
 *  Keeps track of nodes to write data to on upper surface in z
 */
      std::vector<int> upperEdgeNodeNums;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

/**
 * Determines if two nodes have the same configuration space coordinates
 */
      bool sameConfigCoords(int srcIndex, int tarIndex, double dxMin, const Eigen::MatrixXd& nodeList);
  };
}

#endif // LC_SOL_HEAT_FLUX_CALC_H
