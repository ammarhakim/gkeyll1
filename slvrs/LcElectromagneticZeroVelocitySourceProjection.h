/**
 * @file	LcElectromagneticZeroVelocitySourceProjection.h
 *
 * @brief Takes an input A(z) and computes a DG representation
 * for a distribution function with zero mean velocity
 */

#ifndef LC_ELECTROMAGNETIC_ZERO_VELOCITY_SOURCE_PROJECTION_H
#define LC_ELECTROMAGNETIC_ZERO_VELOCITY_SOURCE_PROJECTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using a nodal discontinous
 * Galerkin scheme.
 */
  class ElectromagneticZeroVelocitySourceProjection : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      ElectromagneticZeroVelocitySourceProjection();

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
/** Species mass (kg) */
      double speciesMass;
/** Species temperature (eV) */
      double speciesTemp;
/** Species charge (C) */
      double speciesCharge;
/**
 * Matrix of gaussian quadrature locations.
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing integrals.
 */
      Eigen::MatrixXd gaussOrdinates;
/** Weights for gaussian quadrature points */
      std::vector<double> gaussWeights;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpMatrix;
/** Transpose of interpolation matrix */
      Eigen::MatrixXd interpMatrixTranspose;
/** Inverse of mass matrix */
      Eigen::MatrixXd massMatrixInv;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_ELECTROMAGNETIC_ZERO_VELOCITY_SOURCE_PROJECTION_H
