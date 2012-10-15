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
#include <LcExcept.h>
#include <LcLagrangeTensorBasisCalc.h>
#include <LcMathLib.h>

namespace Lucee
{
  LagrangeTensorBasisCalc::LagrangeTensorBasisCalc(unsigned nd)
    : ndim(nd)
  {
    for (unsigned n=0; n<ndim; ++n)
      numNodes.push_back(0);
  }

  void
  LagrangeTensorBasisCalc::calc(Node_t type, const std::vector<unsigned> nn)
  {
    if (nn.size() != ndim)
    {
      Lucee::Except lce("LagrangeTensorBasisCalc::calc: Must specify number of nodes in each direction.");
      lce << " Only " << nn.size() << " directions specified for NDIM " << ndim << std::endl;
      throw lce;
    }

    for (unsigned i=0; i<ndim; ++i) 
      numNodes[i] = nn[i];

// allocate space to store node locations
    nodeLocs.resize(ndim);
    for (unsigned i=0; i<ndim; ++i)
      nodeLocs[i].loc.resize(numNodes[i]);

// compute location of nodes
    if (type == LagrangeTensorBasisCalc::LOBATTO)
    {
    }
    else if (type == LagrangeTensorBasisCalc::GAUSSIAN)
    {
    }
    else if (type == LagrangeTensorBasisCalc::UNIFORM)
    {
    }
  }
}
