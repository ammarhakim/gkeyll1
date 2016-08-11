/**
 * @file	LcDGDiffusionUpdater1D.h
 *
 * @brief	Updater to evaluate (hyper)diffusion operators using nodal DG
 */

#ifndef LC_DG_DIFFUSION_UPDATER_1D_H
#define LC_DG_DIFFUSION_UPDATER_1D_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to evaluate diffusion operator using a DG scheme. This
 * updater hard-codes a bunch of stencils for various schemes, and can
 * be used to test the efficacy of these schemes.
 */
  class DGDiffusionUpdater1D : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new updater.
 */
      DGDiffusionUpdater1D();

/**
 * Destroy updater.
 */
      ~DGDiffusionUpdater1D();

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
/** Diffusion coefficient */
      double alpha;
/** CFL number */
      double cfl;
/** Flag to indicate if to compute only increments */
      bool onlyIncrement;
/** Polynomial order */
      int polyOrder;
/** Tyep of scheme to use */
      unsigned schemeType;
/** Matrix for current cell */
      Lucee::Matrix<double> iMat;
/** Matrix for left cell */
      Lucee::Matrix<double> lowerMat;
/** Matrix for upper cell */
      Lucee::Matrix<double> upperMat;

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

/**
 * Set stencils for LDG-L scheme 
 *
 * @param dx Grid spacing.  
 */
      void calcLDGLStencil(double dx);

/**
 * Set stencils for LDG-R scheme 
 *
 * @param dx Grid spacing.  
 */
      void calcLDGRStencil(double dx);

/**
 * Set stencils for LDG-S scheme 
 *
 * @param dx Grid spacing.  
 */
      void calcLDGSStencil(double dx);

/**
 * Set stencils for RDG scheme 
 *
 * @param dx Grid spacing.  
 */
      void calcRDGStencil(double dx);
  };
}

#endif // LC_DG_DIFFUSION_UPDATER_1D_H
