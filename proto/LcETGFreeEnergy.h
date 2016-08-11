/**
 * @file	LcETGFreeEnergy.h
 *
 * @brief	Updater to compute the free energy for ETG problem.
 */

#ifndef LC_ETG_FREE_ENERGY_H
#define LC_ETG_FREE_ENERGY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcDynVector.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to integrate nodal DG field over the entire domain.
 */
  class ETGFreeEnergy : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Number of components for coordinate arrays etc. */
      static const unsigned NC4 = Lucee::CDIM<4>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC2 = Lucee::CDIM<2>::N;

/** Create new solver */
      ETGFreeEnergy();

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
/** Temperature of adiabatic species (in eV) */
      double bgAdiabaticTemp;
/** Mass of kinetic species (in kg) */
      double kineticMass;
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
  };
}

#endif // LC_ETG_FREE_ENERGY_H
