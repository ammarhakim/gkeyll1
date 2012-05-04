/**
 * @file	LcNodalPoissonBracketUpdater.h
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

#ifndef LC_NODAL_POISSON_BRACKET_UPDATER_H
#define LC_NODAL_POISSON_BRACKET_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve 2D Poisson bracket operator PDE of the form
 *
 *    da
 *   ---- + {a,b} = 0
 *    dt
 *
 * where {a,b} is the Poisson bracket operator. This updater updates
 * the field 'a' using a first-order forward Euler step with the
 * Poisson bracket computed using a DG scheme.
 */
  class NodalPoissonBracketUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      NodalPoissonBracketUpdater();

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
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Type of interface flux to use */
      unsigned fluxType;

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

/**
 * Struct to hold data for Guassian quadrature.
 */
      struct GaussQuadData
      {
/**
 * Reset object.
 * 
 * @param nord Numer of ordinates in each direction.
 * @param nlocal Total number of local nodes.
 */
          void reset(unsigned nord, unsigned nlocal)
          {
// allocate memory for various matrices
            ords.m = Lucee::Matrix<double>(nord, 3);
            weights.resize(nord);
            interpMat.m = Lucee::Matrix<double>(nord, nlocal);

            for (unsigned dir=0; dir<2; ++dir)
              pDiffMatrix[dir].m = Lucee::Matrix<double>(nord, nlocal);
          }

/** Matrix of ordinates */
          MatrixHolder ords;
/** Vector of weights */
          std::vector<double> weights;
/** Interpolation matrix */
          MatrixHolder interpMat;
/** Differentiation matrices, computing derivatives at quadrature nodes */
          MatrixHolder pDiffMatrix[2];
      };

/** Differentiation matrices */
      MatrixHolder diffMatrix[2];
/** Differentiation matrices */
      MatrixHolder stiffMatrix[2];
/** Liftness matrix on lower edges */
      MatrixHolder lowerLift[2];
/** Liftness matrix on upper edges */
      MatrixHolder upperLift[2];
/** Data for volume quadrature */
      GaussQuadData volQuad;
/** Data for quadrature on each lower face */
      GaussQuadData surfLowerQuad[2];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[2];
/** Gradients of basis functions at quadrature points */
      MatrixHolder mGradPhi[2];
/** Basis fucntions at surface quadrature point */
      MatrixHolder mSurfLowerPhi[2];
/** Basis fucntions at surface quadrature point */
      MatrixHolder mSurfUpperPhi[2];

/**
 * Structure to store node numbers on edges.
 */
      struct EdgeNodeNums
      {
/** Node numbers */
          std::vector<int> nums;
      };

/** Vector to store lower node numbers */
      EdgeNodeNums lowerNodeNums[2];
/** Vector to store upper node numbers */
      EdgeNodeNums upperNodeNums[2];

/**
 * Structure to hold speeds at each node.
 */
      struct NodeSpeed
      {
/** Speeds at each node */
          std::vector<double> s;
      };

/**
 * Calculate speeds in the X and Y directions. The output structure
 * must be pre-allocated.
 * 
 * @param phiK values of potential at nodes.
 * @param speeds On output, speeds in X- and Y-directions.
 */
      void calcSpeeds(std::vector<double>& phiK,
        NodeSpeed speeds[2]);

/**
 * Calculate speeds in the X and Y directions at quadrature
 * points. The output structure must be pre-allocated.
 * 
 * @param phiK values of potential at nodes.
 * @param speeds On output, speeds in X- and Y-directions.
 */
      void calcSpeedsAtQuad(std::vector<double>& phiK,
        NodeSpeed speeds[2]);

/**
 * Return upwind flux based given speed and nodal values.
 *
 * @param u Speed.
 * @param chil Vorticity on left node.
 * @param chir Vorticity on roght node.
 * @rturn upwind flux.
 */
      double getUpwindFlux(double u, double chil, double chir);

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

#endif // LC_NODAL_POISSON_BRACKET_UPDATER_H
