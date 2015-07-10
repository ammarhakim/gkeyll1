/**
 * @file	LcETGInnerProduct.h
 *
 * @brief	Updater to compute the inner product for ETG problem between two distributions.
 */

#ifndef LC_ETG_INNER_PRODUCT_H
#define LC_ETG_INNER_PRODUCT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

namespace Lucee
{
/**
 * Updater to integrate nodal DG field over the entire domain.
 */
  class ETGInnerProduct : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Number of components for coordinate arrays etc. */
      static const unsigned NC4 = Lucee::CDIM<4>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC2 = Lucee::CDIM<2>::N;

/** Create new solver */
      ETGInnerProduct();

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
      Lucee::NodalFiniteElementIfc<4> *nodalBasis4d;
/** Pointer to 2D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Total number of nodes in system */
      int totalNodes;
/** Total number of nodes at a fixed x and y */
      int nodesPerPosition;
/** Reference to function providing write location */
      int fnRef;
/** Temperature of adiabatic species (in eV) */
      double bgAdiabaticTemp;
/** Conversion factor of joules to eV */
      double eV;
/** Density of kinetic species (in 1/m^3) */
      double bgKineticDensity;
/** Mass of kinetic species (in kg) */
      double kineticMass;
/** Stores list of entries to insert into sparse matrix */
      std::vector<Eigen::Triplet<double> > tripletList;
/** Name of output file */
      std::string filename;
/**
 * Struct to hold data for Guassian quadrature.
 */
  struct GaussQuadData
      {
/**
 * Reset object.
 * 
 * @param numQuad Numer of quadrature points.
 * @param nlocal Total number of local nodes.
 */
          void reset(int numQuad, int nlocal)
          {
            // allocate memory for various matrices
            weights = Eigen::VectorXd(numQuad);
            interpMat = Eigen::MatrixXd(numQuad, nlocal);
          }

          /** Vector of weights */
          Eigen::VectorXd weights;
          /** Interpolation matrix */
          Eigen::MatrixXd interpMat;
      };

      GaussQuadData volQuad2d;
      GaussQuadData volQuad4d;

/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
/**
 * Use to figure out what two nodes are nonzero
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        std::vector<double>& res);
  };
}

#endif // LC_ETG_INNER_PRODUCT_H
