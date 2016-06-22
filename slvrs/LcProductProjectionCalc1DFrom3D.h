/**
 * @file	LcProductProjectionCalc1DFrom3D.h
 *
 * @brief	Projects the product of g*f onto a 1d field, where g and f are 3d fields and g is a 1d field
 */

#ifndef LC_PRODUCT_PROJECTION_CALC_1D_FROM_3D_H
#define LC_PRODUCT_PROJECTION_CALC_1D_FROM_3D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Applies particle refection BCs to distribution function
 */
  class ProductProjectionCalc1DFrom3D : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new projection updater */
      ProductProjectionCalc1DFrom3D();
/**  Destructor */
      virtual ~ProductProjectionCalc1DFrom3D();

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
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to configuration space basis functions */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis1d;
/** Factor to multiply all results by (like 2*pi*B/m to account v_perp -> mu integration */
      double scaleFactor;
/** Space to store moment data for on a processor */
      Lucee::Field<1, double> *resultLocal;
/**
 * Interpolation matrix for 3d quadrature
 */
      Eigen::MatrixXd interpMatrix3d;
/**
 * Interpolation matrix for 1d quadrature
 */
      Eigen::MatrixXd interpMatrix1d;
/**
 * Inverse of 1d mass matrix
 */
      Eigen::MatrixXd massMatrixInv1d;
/**
 * Gaussian quadrature weights for 3d integration
 */
      std::vector<double> gaussWeights3d;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_PRODUCT_PROJECTION_CALC_1D_FROM_3D_H
