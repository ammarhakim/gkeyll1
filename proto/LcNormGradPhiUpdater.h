/**
 * @file	LcNormGradPhiUpdater.h
 *
 * @brief	Updater to compute |grad.p|^2 integrated over the domain.
 */

#ifndef LC_NORM_GRAD_PHI_UPDATER_H
#define LC_NORM_GRAD_PHI_UPDATER_H

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
 * Updater to |grad.p|^2 integrated over the domain.
 */
  template <unsigned NDIM>
  class NormGradPhiUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Ctor */
      NormGradPhiUpdater();

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
      MatrixHolder diffMatrix[NDIM];
/** Interpolation matrix */
      MatrixHolder interpMat;
/** Weights for quadrature */
      std::vector<double> weights;
/** Ordinates for quadrature */
      MatrixHolder ordinates;
/** Differentiation matrices, computing derivatives at quadrature nodes */
      MatrixHolder pDiffMatrix[NDIM];

/**
 * Compute norm of gradient of field at each node.
 *
 * @param phiK potential at nodes
 * @param normGradPhi Norm of gradient phi.
 */
      void calcNormGrad(std::vector<double>& phiK,
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
        const double* vec, double v, double* out);
  };
}

#endif // LC_NORM_GRAD_PHI_UPDATER_H
