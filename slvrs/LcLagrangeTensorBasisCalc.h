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
/**
 * Create new tensor basis calculator object.
 *
 * @param n Dimension of space.
 */
      LagrangeTensorBasisCalc(unsigned n);

    private:
/** Dimension of basis functions */
      unsigned ndim;
/** Number of basis functions in each direction */
      std::vector<unsigned> numBasis;
  };
}

#endif // LC_LAGRANGE_TENSOR_BASIS_CALC_H
