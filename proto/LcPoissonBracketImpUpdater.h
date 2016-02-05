/**
 * @file	LcPoissonBracketImpUpdater.h
 *
 * @brief	Updater to solver Poisson bracket operator PDEs.
 * This version is used to test speed improvements to the existing updater.
 */

#ifndef LC_POISSON_BRACKET_IMP_UPDATER_H
#define LC_POISSON_BRACKET_IMP_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcPoissonBracketEquation.h>

// eigen includes
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
  class PoissonBracketImpUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC = Lucee::CDIM<NDIM>::N;

/** Create new nodal DG solver */
      PoissonBracketImpUpdater();

/** Dtor */
      virtual ~PoissonBracketImpUpdater();      

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
/** Equation to solve */
      Lucee::PoissonBracketEquation *equation;
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
          void reset(int numQuad, int nlocal, int numCoords)
          {
            // allocate memory for various matrices
            coordMat = Eigen::MatrixXd(numQuad, numCoords);
            weights = Eigen::VectorXd(numQuad);
            interpMat = Eigen::MatrixXd(numQuad, nlocal);
            nodeNums.resize(nlocal);
          }

/** Matrix of ordinates (could remove this) */
          Eigen::MatrixXd coordMat;
/** Vector of weights */
          Eigen::VectorXd weights;
/** Interpolation matrix */
          Eigen::MatrixXd interpMat;
/** Vector of node numbers (only needed for surface structures) */
          std::vector<int> nodeNums;
      };
/** Data for volume quadrature */
      GaussQuadData volQuad;
/** Data for quadrature on each lower face */
      GaussQuadData surfLowerQuad[NDIM];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[NDIM];
      std::vector<Eigen::MatrixXd> gradMatrices;
/** Inverse of mass matrix */
      Eigen::MatrixXd massMatrixInv;
/** Stored matrices for stupid testing */
      std::vector<std::vector<Eigen::MatrixXd> > bigStoredUpperSurfMatrices;
      std::vector<std::vector<Eigen::MatrixXd> > bigStoredLowerSurfMatrices;
      std::vector<std::vector<Eigen::MatrixXd> > bigStoredVolMatrices;
/** Flag to indicate if only increments should be computed */
      bool onlyIncrement;
/** Flag to indicate if a Jacobian factor is supplied */
      bool hasJacobian;
/** Field to store optional Jacobian field (SHOULD REALLY BE SPATIAL DIM SIZE) */
      Lucee::Field<NDIM, double> *jacobianField;
/** Impional vector storing directions to update */
      std::vector<int> updateDirs;
/** Impional vector keeping track of whether to apply zero flux BCs in each direction */
      std::vector<bool> zeroFluxFlags;
 /**
 * CompsurfLowerQuad[dir].interpMat*rightDatasurfLowerQuad[dir].interpMat*rightDataute numerical flux
 * @param alphaDotN: characteristic velocities at quad points (dot n)
 * @param leftValsAtQuad: left cell values at quad points
 * @param rightValsAtQuad: right cell values at quad points
 * @param numericalFluxAtQuad: resulting numerical flux calculation at quad points
 */
      void computeNumericalFlux(const Eigen::VectorXd& alphaDotN, const Eigen::VectorXd& leftValsAtQuad,
          const Eigen::VectorXd& rightValsAtQuad, Eigen::VectorXd& numericalFluxAtQuad);
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

// Stuff for timing studies
      double totalVolTime, totalSurfTime, jacAtQuad;
      double computeVolAlphaAtQuadNodesTime, computeSurfAlphaAtQuadNodesTime;
      double vol_loop1, vol_loop2, vol_loop3, vol_loop4;
      double surf_loop1, surf_loop2, surf_loop3, surf_loop4;
  };
}

#endif // LC_POISSON_BRACKET_IMP_UPDATER_H
