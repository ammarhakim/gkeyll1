/**
 * @file	LcETGAdjointSource.h
 *
 * @brief	Updater to implement adjoint gyrokinetic equation
 * for the ETG problem
 */

#ifndef LC_ETG_ADJOINT_SOURCE_H
#define LC_ETG_ADJOINT_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/LU>
#include <Eigen/Dense>

namespace Lucee
{
  template <unsigned NDIM>
  class ETGAdjointSource : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC = Lucee::CDIM<NDIM>::N;

/** Create new nodal DG solver */
      ETGAdjointSource();

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

      double kineticMass;
      double eV;
      double tauOverZi;
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
            pDiffMatrix.resize(nlocal);
            interpMat = Eigen::MatrixXd(numQuad, nlocal);

            for (int i = 0; i < nlocal; i++)
              pDiffMatrix[i] = Eigen::MatrixXd(NDIM, numQuad);
          }

/** Matrix of ordinates (could remove this) */
          Eigen::MatrixXd coordMat;
/** Vector of weights */
          Eigen::VectorXd weights;
/** Interpolation matrix */
          Eigen::MatrixXd interpMat;
/** Differentiation matrices, computing derivative of each basis function at quad points */
          std::vector<Eigen::MatrixXd> pDiffMatrix;
      };
/** Data for volume quadrature */
      GaussQuadData volQuad;
/** Data for quadrature on each lower face */
      GaussQuadData surfLowerQuad[NDIM];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[NDIM];
/** Inverse of mass matrix */
      Eigen::MatrixXd massMatrixInv;
/** Stored matrices for stupid testing */
      std::vector<Eigen::MatrixXd> bigStoredVolMatrices;
/** Flag to indicate if only increments should be computed */
      bool onlyIncrement;
/** Flag to indicate if a Jacobian factor is supplied */
      bool hasJacobian;
/** Field to store optional Jacobian field (SHOULD REALLY BE SPATIAL DIM SIZE) */
      Lucee::Field<NDIM, double> *jacobianField;
/** Optional vector storing directions to update */
      std::vector<int> updateDirs;
/** Optional vector keeping track of whether to apply zero flux BCs in each direction */
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
  };
}

#endif // LC_ETG_ADJOINT_SOURCE_H
