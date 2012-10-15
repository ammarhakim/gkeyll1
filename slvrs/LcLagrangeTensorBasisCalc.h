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

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to calculate data needed for Lagrange tensor basis functions.
 */
  class LagrangeTensorBasisCalc
  {
    public:
/** Enum for location of nodes */
      enum Node_t { LOBATTO, GAUSSIAN, UNIFORM };

/**
 * Create new tensor basis calculator object.
 *
 * @param n Dimension of space.
 */
      LagrangeTensorBasisCalc(unsigned n);

/**
 *  Calculate the data needed in constructing the basis
 *  function. Calling this function again will redo the calculation.
 *
 * @param type Type of nodal layout.
 * @param numNodes Number of nodes in each direction. Should have exactly ndim elements.
 */
      void calc(Node_t type, const std::vector<unsigned> numNodes);

    private:
/** Dimension of basis functions */
      unsigned ndim;
/** Number of nodes in each direction */
      std::vector<unsigned> numNodes;

/** Structure to store location of nodes */
      struct NodeLoc
      {
/** Location of nodes */
          std::vector<double> loc;
      };

/** Location of nodes in each direction */
      std::vector<NodeLoc> nodeLocs;
  };
}

#endif // LC_LAGRANGE_TENSOR_BASIS_CALC_H
