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
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcMatrix.h>
#include <LcRowMajorIndexer.h>
#include <LcRowMajorSequencer.h>

namespace Lucee
{
  template <unsigned NDIM>
  LagrangeTensorBasisCalc<NDIM>::LagrangeTensorBasisCalc()
  {
    for (unsigned n=0; n<NDIM; ++n)
      numNodes.push_back(0);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::calc(Node_t type, const std::vector<unsigned> nn)
  {
    if (nn.size() != NDIM)
    {
      Lucee::Except lce("LagrangeTensorBasisCalc::calc: Must specify number of nodes in each direction.");
      lce << " Only " << nn.size() << " directions specified for NDIM " << NDIM << std::endl;
      throw lce;
    }

    for (unsigned i=0; i<NDIM; ++i) 
      numNodes[i] = nn[i];

// allocate space to store node locations
    nodeLocs.resize(NDIM);
    for (unsigned i=0; i<NDIM; ++i)
      nodeLocs[i].loc.resize(numNodes[i]);

// compute location of nodes
    if (type == Lucee::LOBATTO)
      createLobattoNodes();
    else if (type == Lucee::GAUSSIAN)
      createGaussianNodes();
    else if (type == Lucee::UNIFORM)
      createUniformNodes();

// total number of nodes in element
    unsigned totalNodes = 1;
    for (unsigned d=0; d<NDIM; ++d)
      totalNodes *= numNodes[d];

    int lower[NDIM], upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
    {
      lower[d] = 0;
      upper[d] = numNodes[d];
    }
    Lucee::Region<NDIM, int> rgn(lower, upper);
// create sequencer for looping over nodes (why RowMajor and not
// ColumnMajor? It actually does not matter, as long as all other
// calculations are done consistently.).
    Lucee::RowMajorSequencer<NDIM> nodeSeq(rgn), polySeq(rgn);
// indexer to map to location in matrix
    Lucee::RowMajorIndexer<NDIM> idx(rgn);

    double xn[NDIM];
// create matrix of linear system coefficients: each node contributes
// one row. In that row each entry is the product of Legendre
// polynomials evaluated at that node.
    Lucee::Matrix<double> coeffMat(totalNodes, totalNodes);

    int nodeIdx[NDIM], polyIdx[NDIM];
    while (nodeSeq.step())
    {
      nodeSeq.fillWithIndex(nodeIdx);
      int row = idx.getIndex(nodeIdx); // row of coeff matrix

// compute nodal coordinates: this is not strictly needed, but makes
// following code less confusing.
      for (unsigned d=0; d<NDIM; ++d)
        xn[d] = nodeLocs[d].loc[nodeIdx[d]];

// each entry in this row is product of Legendre polynomials evaluate
// with nodal coordinate components
      polySeq.reset();
      while (polySeq.step())
      {
        polySeq.fillWithIndex(polyIdx);
        int col = idx.getIndex(polyIdx);

        double entry = 1.0;
        for (unsigned d=0; d<NDIM; ++d)
          entry *= Lucee::legendrePoly(polyIdx[d], xn[d]);
// insert this into appropriate location in matrix
        coeffMat(row, col) = entry;
      }
    }

// contruct RHS matrix, which is just a unit matrix as basis functions
// are such that they evaluate to 1 at one node and are zero at other
// nodes.
    expandCoeff = Lucee::Matrix<double>(totalNodes, totalNodes);
    expandCoeff = 0.0;
    for (unsigned i=0; i<totalNodes; ++i)
      expandCoeff(i,i) = 1.0;

// invert system to get expansion coefficients
    Lucee::solve(coeffMat, expandCoeff);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::createLobattoNodes()
  {
    throw Lucee::Except("LagrangeTensorBasisCalc::createLobattoNodes: Not implemented!");
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::createGaussianNodes()
  {
    for (unsigned d=0; d<NDIM; ++d)
    {
      Lucee::Vector<double> x(numNodes[d]), w(numNodes[d]);
      Lucee::gauleg(numNodes[d], -1, 1, x, w); // nodes at Gaussian quadrature points
      for (unsigned i=0; i<numNodes[d]; ++i)
        nodeLocs[d].loc[i] = x[i];
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::createUniformNodes()
  {
    for (unsigned d=0; d<NDIM; ++d)
    {
      double dx = 2.0/(numNodes[d]-1);
      nodeLocs[d].loc[0] = -1.0; // first node on left edge
      for (unsigned i=1; i<numNodes[d]; ++i)
        nodeLocs[d].loc[i] += nodeLocs[d].loc[i]+dx; // nodes are evenly spaced
    }
  }

// instantiations
  template class LagrangeTensorBasisCalc<1>;
  template class LagrangeTensorBasisCalc<2>;
  template class LagrangeTensorBasisCalc<3>;
  template class LagrangeTensorBasisCalc<4>;
  template class LagrangeTensorBasisCalc<5>;
}
