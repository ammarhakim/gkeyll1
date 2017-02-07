/**
 * @file	LcLenardBernsteinDiffAlternate3D2VUpdater.h
 *
 * @brief	Updater to evaluate the diffusion term in the L-B collision operator for 3D2V problems.
 * This updater is different from LcLenardBernsteinDiff5DUpdater in its use of a 2D recovery method
 * to evaluate diffusion terms
 */

#ifndef LC_LENARD_BERNSTEIN_DIFF_ALTERNATE_3D_2V_UPDATER_H
#define LC_LENARD_BERNSTEIN_DIFF_ALTERNATE_3D_2V_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcStructGridField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to evaluate the diffusion term in the Lenard-Bernstein
 * collision operator.
 */
  class LenardBernsteinDiffAlternate3D2VUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC2 = Lucee::CDIM<2>::N;

/**
 * Create new updater.
 */
      LenardBernsteinDiffAlternate3D2VUpdater();

/**
 * Destroy updater.
 */
      ~LenardBernsteinDiffAlternate3D2VUpdater();

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
/** CFL number */
      double cfl;
/** Flag to indicate if to compute only increments */
      bool onlyIncrement;
/** Reference to function specifying alpha */
      int fnRef;
/** mass of species */
      double speciesMass;
/** Polynomial order of basis functions */
      int polyOrder;
/** Store a vector of 1d quadrature points on reference element */
      std::vector<double> gaussPoints1d;
/** Store a vector of recovery matrices for each surface in v_parallel */
      std::vector<Eigen::MatrixXd> recoveryMatricesVPara;
/** Store a vector of recovery matrices for each surface in mu */
      std::vector<Eigen::MatrixXd> recoveryMatricesMu;
/** Store interpolation matrix for recovery procedure in v_parallel */
      std::vector<Eigen::MatrixXd> recoveryInterpMatVPara;
/** Store interpolation matrix for recovery procedure in mu*/
      std::vector<Eigen::MatrixXd> recoveryInterpMatMu;
/** Pointer to 5d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<5> *nodalBasis5d;
/** Pointer to 3d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 2d nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Keeps track of the offsets needed to get all nodes that share the same config. space location */
      std::vector<int> nodalStencil;
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
      GaussQuadData surfLowerQuad[2];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[2];
/** List of matrices on each lower face */
      std::vector<Eigen::MatrixXd > lowerMat;
/** List of matrices on each upper face */
      std::vector<Eigen::MatrixXd > upperMat;
/** Precomputed matrices accounting for how mu is a radial coordinate */
      Eigen::MatrixXd lowerMatTimesMu;
      Eigen::MatrixXd upperMatTimesMu;
/** Weights for volume gaussian quadrature points */
      std::vector<double> gaussVolWeights;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to volume
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpVolMatrix;
/**
 * Matrix of surface gaussian quadrature locations on bottom face..
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing surface integrals.
 */
      Eigen::MatrixXd gaussSurfOrdinates;
/** Weights for surface gaussian quadrature points */
      std::vector<double> gaussSurfWeights;
/** 
 * Interpolation matrix for bringing data that lives on surface in diffDir
 * to quadrature points on same surface. It is constructed from the bottom surf
 * basis function evaluations but should work for upper surface as well.
 * The quadrature locations for these nodes are also in gaussSurfOrdinates and
 * gaussSurfWeights.
 */
      Eigen::MatrixXd surfNodeInterpMatrix;

      std::vector<Eigen::MatrixXd> surfLowerUpdateMatrix;
      std::vector<Eigen::MatrixXd> surfUpperUpdateMatrix;
      std::vector<Eigen::MatrixXd> surfLowerUpdateMatrixDeriv;
      std::vector<Eigen::MatrixXd> surfUpperUpdateMatrixDeriv;
      std::vector<Eigen::MatrixXd> volumeUpdateMatrix;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
/**
 * Evaluate function at specified location and fill output array with
 * result.
 *
 * @param L Lua state object to use.
 * @param tm Time to evaluate function at.
 * @param res On output, result of evaluating function.
 */
      void evaluateFunction(Lucee::LuaState& L, double tm, 
        std::vector<double>& res);
/**
 * Determines if two nodes have the same configuration space coordinates
 */
      bool sameConfigCoords(int srcIndex, int tarIndex, double dxMin, const Eigen::MatrixXd& nodeList);
  };
}

#endif // LC_LENARD_BERNSTEIN_DIFF_ALTERNATE_3D_2V_UPDATER_H
