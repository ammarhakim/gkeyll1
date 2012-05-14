/**
 * @file	LcEnergyFromStreamFunctionUpdater.h
 *
 * @brief	Updater to compute energy from streamfunction.
 */

#ifndef LC_ENERGY_FROM_STREAM_FUNCTION_UPDATER_H
#define LC_ENERGY_FROM_STREAM_FUNCTION_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDynVector.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to compute total energy in domain given the streamfunction.
 */
  class EnergyFromStreamFunctionUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      EnergyFromStreamFunctionUpdater();

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

/**
 * Matrix holder: this class is needed as the Matrix class does not
 * have a default constructor.
 */
      struct MatrixHolder
      {
/** Ctor */
          MatrixHolder() : m(1, 1) {}
/** Differentiation matrix */
          Lucee::Matrix<double> m;
      };

/** Differentiation matrices */
      MatrixHolder diffMatrix[2];
/** Interpolation matrix */
      MatrixHolder interpMat;
/** Weights for quadrature */
      std::vector<double> weights;
/** Ordinates for quadrature */
      MatrixHolder ordinates;
/** Differentiation matrices, computing derivatives at quadrature nodes */
      MatrixHolder pDiffMatrix[2];

/**
 * Compute norm of gradient of field at each node.
 *
 * @param phiK potential at nodes
 * @param normGradPhi Norm of gradient phi.
 */
      void calcNormGrad(std::vector<double>& phiK,
        std::vector<double>& normGradPhi);

/**
 * Compute norm of gradient of field at each node.
 *
 * @param phiK potential at nodes
 * @param normGradPhi Norm of gradient phi.
 */
      void calcNormGrad_1(std::vector<double>& phiK,
        std::vector<double>& normGradPhi);

/**
 * Compute matrix-vector multiply. Output vector must be
 * pre-allocated. Note that the computation performed is
 *
 * out = m*mat*vec + v*out
 *
 * @param m Factor for accumulation.
 * @param mat Matrix for multiplication.
 * @param vec Vector for multiplication.
 * @param v Factor for accumulation.
 * @param out On output, holds the product.
 */
      void matVec(double m, const Lucee::Matrix<double>& mat,
        const std::vector<double>& vec, double v, double* out);
  };
}

#endif // LC_ENERGY_FROM_STREAM_FUNCTION_UPDATER_H
