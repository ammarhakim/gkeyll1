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

    private:
/** Number of nodes in each direction */
      unsigned numNodes[NDIM];
/** Total number of nodes */
      unsigned totalNodes;

/** Structure to store location of nodes */
      struct NodeLoc
      {
/** Location of nodes */
          std::vector<double> loc;
      };

/** Location of nodes in each direction */
      NodeLoc nodeLocs[NDIM];

/** Matrix of expansion coefficients */
      Lucee::Matrix<double> expandCoeff;

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
  };
}

#endif // LC_LAGRANGE_TENSOR_BASIS_CALC_H
