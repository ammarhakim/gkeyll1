/**
 * @file	LcPoissonBracketUpdater.h
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 */

#ifndef LC_POISSON_BRACKET_UPDATER_H
#define LC_POISSON_BRACKET_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcHyperEquation.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve general dimension Poisson bracket operator PDE of the form
 *
 *    da
 *   ---- + {a,b} = 0
 *    dt
 *
 * where {a,b} is the Poisson bracket operator. This updater updates
 * the field 'a' using a first-order forward Euler step with the
 * Poisson bracket computed using a DG scheme.
 */
  template <unsigned NDIM>
  class PoissonBracketUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      PoissonBracketUpdater();

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
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Type of interface flux to use */
      unsigned fluxType;
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
            coordMat = Eigen::MatrixXd(numQuad, NDIM);
            weights = Eigen::VectorXd(numQuad);
            pDiffMatrix.resize(nlocal);
            interpMat = Eigen::MatrixXd(numQuad, nlocal);

            for (int dir = 0; dir < NDIM; dir++)
              pDiffMatrix[dir] = Eigen::MatrixXd(NDIM, numQuad);
          }

/** Matrix of ordinates */
          Eigen::MatrixXd coordMat;
/** Vector of weights */
          Eigen::VectorXd weights;
/** Interpolation matrix */
          Eigen::MatrixXd interpMat;
/** Differentiation matrices, computing derivative of each basis function at quad points */
          std::vector<Eigen::MatrixXd> pDiffMatrix;
      };

/** Differentiation matrices */
      std::vector<Eigen::MatrixXd> diffMatrix;
/** Differentiation matrices */
      std::vector<Eigen::MatrixXd> gradStiffnessMatrix;
/** Data for volume quadrature */
      GaussQuadData volQuad;
/** Data for quadrature on each lower face */
      GaussQuadData surfLowerQuad[NDIM];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[NDIM];
/** Gradients of basis functions at quadrature points */
      std::vector<Eigen::MatrixXd> mGradPhi;
/** Basis fucntions at surface quadrature point */
      std::vector<Eigen::MatrixXd> mSurfLowerPhi;
/** Basis fucntions at surface quadrature point */
      std::vector<Eigen::MatrixXd> mSurfUpperPhi;
/** Inverse of mass matrix */
      Eigen::MatrixXd massMatrixInv;
/** Flag to indicate if only increments should be computed */
      bool onlyIncrement;
/** Vector to store lower node numbers */
      std::vector<std::vector<int> > lowerNodeNums;
/** Vector to store upper node numbers */
      std::vector<std::vector<int> > upperNodeNums;

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
        NodeSpeed speeds[NDIM]);

/**
 * Calculate speeds in the X and Y directions at quadrature
 * points. The output structure must be pre-allocated.
 * 
 * @param phiK values of potential at nodes.
 * @param speeds On output, speeds in X- and Y-directions.
 */
      void calcSpeedsAtQuad(std::vector<double>& phiK,
        NodeSpeed speeds[NDIM]);

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
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_POISSON_BRACKET_UPDATER_H
