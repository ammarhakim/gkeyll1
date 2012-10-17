/**
 * @file	LcLagrangeTensorBasisCalc.h
 *
 * @brief	Class to calculate data needed for Lagrange tensor basis functions.
 */
#ifndef LC_LAGRANGE_TENSOR_BASIS_CALC_H
#define LC_LAGRANGE_TENSOR_BASIS_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// gkeyll includes
#include <LcMatrix.h>
#include <LcRowMajorIndexer.h>

// std includes
#include <vector>

namespace Lucee
{
/** Enum for location of nodes */
  enum Node_t { LOBATTO, GAUSSIAN, UNIFORM };

/**
 * Class to calculate data needed for Lagrange tensor basis functions.
 */
  template <unsigned NDIM>
  class LagrangeTensorBasisCalc
  {
    public:

/**
 * Create new tensor basis calculator object.
 */
      LagrangeTensorBasisCalc();

/**
 * Get total number of nodes in element.
 *
 * @return Total number of nodes in element.
 */
      unsigned getNumNodes() const 
      { return totalNodes; }

/**
 *  Calculate the data needed in constructing the basis
 *  function. Calling this function again will redo the calculation.
 *
 * @param type Type of nodal layout.
 * @param numNodes Number of nodes in each direction. Should have exactly ndim elements.
 */
      void calc(Node_t type, const unsigned numNodes[NDIM]);

/**
 * Get node location in specified direction.
 *
 * @param dir Direction for node locations.
 * @return vector of node locations.
 */
      std::vector<double> getNodeLoc(unsigned dir) const;

/**
 * Get coordinate of specified node.
 *
 * @param nIdx Node index (0-based)
 * @param xn On output, coordinate of node.
 */
      void fillWithNodeCoordinate(unsigned nIdx, double xn[NDIM]) const;

/**
 * Fetch coefficient matrix.
 *
 * @param coeff On output, this contains the coefficient matrix. Should be pre-allocated.
 */
      void getCoeffMat(Lucee::Matrix<double>& coeff) const;

/**
 * Fetch mass-matrix.
 *
 * @param mMatrix On output, this contains the mass matrix. Should be pre-allocated.
 */
      void getMassMatrix(Lucee::Matrix<double>& mMatrix) const;

/**
 * Evaluate specified basis function at location. This method is very
 * slow, and so should not be called inside an inner loop. It is
 * provided to allow initialization of data needed in updaters, so its
 * cost will be amortized.
 *
 * @param bIdx Basis function index (0-based).
 * @param xc Coorinates in element.
 * @return value of basis function at location.
 */
      double evalBasis(unsigned bIdx, double xc[NDIM]) const;

    private:
/** Number of nodes in each direction */
      unsigned numNodes[NDIM];
/** Total number of nodes */
      unsigned totalNodes;
/** Region representing node shape */
      typename Lucee::Region<NDIM, int> nodeRgn;

/** Structure to store location of nodes */
      struct NodeLoc
      {
/** Location of nodes */
          std::vector<double> loc;
      };

/** Location of nodes in each direction */
      NodeLoc nodeLocs[NDIM];

/** Structure to store coordinates of each node */
      struct NodeCoord
      {
/** Node coordinate */
          double x[NDIM];
      };

/** Coordinates of each node (this seems redundant, but is stored
 * anyway to make life easier) */
      std::vector<NodeCoord> nodeCoords;

/** Matrix of expansion coefficients */
      Lucee::Matrix<double> expandCoeff;

/** Mass matrix for element */
      Lucee::Matrix<double> massMatrix;

/**
 * Create nodes located at Lobatto quadrature points.
 */
      void createLobattoNodes();

/**
 * Create nodes located at Gaussian quadrature points.
 */
      void createGaussianNodes();

/**
 * Create nodes with uniform spacing.
 */
      void createUniformNodes();

/**
 * Compute mass matrix.
 */
      void calcMassMatrix();
  };
}

#endif // LC_LAGRANGE_TENSOR_BASIS_CALC_H
