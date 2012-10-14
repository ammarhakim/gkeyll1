/**
 * @file	LcLagrangeTensorBasisCalc.cpp
 *
 * @brief	Class to calculate data needed for Lagrange tensor basis functions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// gkeyll includes
#include <LcLagrangeTensorBasisCalc.h>

namespace Lucee
{
  LagrangeTensorBasisCalc(unsigned nd)
    : ndim(nd)
  {
// allocate space to store number of basis in each direction
    for (unsigned n=0; n<ndim; ++n)
      numBasis.push_back(0);
  }
}
