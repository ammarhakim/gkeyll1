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
    : nodeLayout(Lucee::UNIFORM)
  {
// by default, a single node
    unsigned nn[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
      nn[d] = 1;
    calc(Lucee::UNIFORM, nn);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::calc(Node_t type, const unsigned nn[NDIM])
  {
// copy over basic data
    nodeLayout = type;
    for (unsigned i=0; i<NDIM; ++i) 
      numNodes[i] = nn[i];

// allocate space to store node locations
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
    totalNodes = 1;
    for (unsigned d=0; d<NDIM; ++d)
      totalNodes *= numNodes[d];

    int lower[NDIM], upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
    {
      lower[d] = 0;
      upper[d] = numNodes[d];
    }
    Lucee::Region<NDIM, int> rgn(lower, upper);
    nodeRgn = rgn;
// create sequencer for looping over nodes (why RowMajor and not
// ColumnMajor? It actually does not matter, as long as all other
// calculations are done consistently. However, it does matter in
// interpreting the output from an updater, as the data is arranged
// with the same layout as used here.).
    Lucee::RowMajorSequencer<NDIM> nodeSeq(rgn), polySeq(rgn);
// indexer to map to location in matrix
    Lucee::RowMajorIndexer<NDIM> idx(rgn);

    int nodeIdx[NDIM], polyIdx[NDIM];

// store nodal coordinates
    nodeCoords.resize(totalNodes);
    while (nodeSeq.step())
    {
      nodeSeq.fillWithIndex(nodeIdx);
      int nn = idx.getIndex(nodeIdx); // node number

      for (unsigned d=0; d<NDIM; ++d)
        nodeCoords[nn].x[d] = nodeLocs[d].loc[nodeIdx[d]];
    }

    double xn[NDIM];
// create matrix of linear system coefficients: each node contributes
// one row. In that row each entry is the product of Legendre
// polynomials evaluated at that node.
    Lucee::Matrix<double> coeffMat(totalNodes, totalNodes);

    nodeSeq.reset();
    while (nodeSeq.step())
    {
      nodeSeq.fillWithIndex(nodeIdx);
      int row = idx.getIndex(nodeIdx); // row of coeff matrix

// each entry in this row is product of Legendre polynomials evaluate
// with nodal coordinate components
      polySeq.reset();
      while (polySeq.step())
      {
        polySeq.fillWithIndex(polyIdx);
        int col = idx.getIndex(polyIdx);

// compute product of Legendre polynomials
        double entry = 1.0;
        for (unsigned d=0; d<NDIM; ++d)
          entry *= Lucee::legendrePoly(polyIdx[d], nodeCoords[row].x[d]);

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

// compute various matrices and quadrature data
    calcMassMatrix();
    for (unsigned d=0; d<NDIM; ++d)
      calcGradStiff(d);

// clear out existing data in node lists
    exclNodes.clear();
    for (unsigned d=0; d<NDIM; ++d)
    {
      lowerNodes[d].nodes.clear();
      upperNodes[d].nodes.clear();
    }

// compute node lists
    if (type == Lucee::GAUSSIAN)
    {
// there are no surface nodes for Gaussian elements
      numExclNodes = nodeRgn.getVolume();

// compute exclusively owned nodes
      Lucee::RowMajorSequencer<NDIM> seq(nodeRgn);
      Lucee::RowMajorIndexer<NDIM> indexer(nodeRgn);

      int idx[NDIM];      
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        exclNodes.push_back(indexer.getIndex(idx));
      }
    }
    else
    { // LOBATTO and UNIFORM
      for (unsigned d=0; d<NDIM; ++d)
      {
        lower[d] = 0;
        upper[d] = numNodes[d]-1;
      }
      numExclNodes = Lucee::Region<NDIM, int>(lower, upper).getVolume();

// compute exclusively owned nodes
      Lucee::RowMajorSequencer<NDIM> seq(nodeRgn);
      Lucee::RowMajorIndexer<NDIM> indexer(nodeRgn);

      int idx[NDIM];
      while (seq.step())
      {
        seq.fillWithIndex(idx);
// ensure node is not on any of the upper faces
        bool onUpper = false;
        for (unsigned d=0; d<NDIM; ++d)
        {
          if ((nodeRgn.getUpper(d)-1) == idx[d])
          {
            onUpper = true;
            break;
          }
        }

        if (!onUpper)
          exclNodes.push_back(indexer.getIndex(idx));
      }

// compute nodes on lower faces
      for (unsigned d=0; d<NDIM; ++d)
      {
// create box only containing face nodes
        Lucee::Region<NDIM, int> faceRgn = nodeRgn.resetBounds(d, nodeRgn.getLower(d), nodeRgn.getLower(d)+1);
       
        seq.reset();
        while (seq.step())
        {
          seq.fillWithIndex(idx);
// add in appropriate list if node is on face
          if (faceRgn.isInside(idx))
            lowerNodes[d].nodes.push_back(indexer.getIndex(idx));
        }
      }

// compute nodes on upper faces
      for (unsigned d=0; d<NDIM; ++d)
      {
// create box only containing face nodes
        Lucee::Region<NDIM, int> faceRgn = nodeRgn.resetBounds(d, nodeRgn.getUpper(d)-1, nodeRgn.getUpper(d));
       
        seq.reset();
        while (seq.step())
        {
          seq.fillWithIndex(idx);
// add in appropriate list if node is on face
          if (faceRgn.isInside(idx))
            upperNodes[d].nodes.push_back(indexer.getIndex(idx));
        }
      }
    }
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
// NOTE: It is perhaps important to ensure exact symmetry of node
// locations (if applicable) as far as possible. Even small errors
// might lead to subtle problems down-stream in updaters that use
// these values. For now, this is a naive implementation, but
// something to keep in mind for the future.

    for (unsigned d=0; d<NDIM; ++d)
    {
      if (numNodes[d] == 1)
        nodeLocs[d].loc[0] = 0.0;
      else
      {
        double dx = 2.0/(numNodes[d]-1);
        nodeLocs[d].loc[0] = -1.0; // first node on left edge
        for (unsigned i=1; i<numNodes[d]-1; ++i)
          nodeLocs[d].loc[i] += nodeLocs[d].loc[i-1]+dx; // nodes are evenly spaced
        nodeLocs[d].loc[numNodes[d]-1] = 1.0; // last node on right edge
      }
    }
  }

  template <unsigned NDIM>
  double
  LagrangeTensorBasisCalc<NDIM>::evalBasis(unsigned bIdx, double xc[NDIM]) const
  {
    Lucee::RowMajorSequencer<NDIM> seq(nodeRgn);
    Lucee::RowMajorIndexer<NDIM> idx(nodeRgn);
    int nodeIdx[NDIM];

    double v = 0;
    while (seq.step())
    {
      seq.fillWithIndex(nodeIdx);
      int nn = idx.getIndex(nodeIdx); // node number

      double pt = expandCoeff(nn, bIdx); // coefficient of expansion
      for (unsigned d=0; d<NDIM; ++d)
        pt *= Lucee::legendrePoly(nodeIdx[d], xc[d]);
      v += pt; // increment sum
    }
    return v;
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::calcMassMatrix()
  {
    Lucee::RowMajorSequencer<NDIM> seq(nodeRgn);
    Lucee::RowMajorIndexer<NDIM> idx(nodeRgn);
    int nodeIdx[NDIM];

// space for mass-matrix
    massMatrix = Lucee::Matrix<double>(totalNodes, totalNodes);

// compute each entry in matrix
    for (unsigned k=0; k<totalNodes; ++k)
    {
      for (unsigned m=0; m<totalNodes; ++m)
      {
        double entry = 0.0;
// loop over basis function expansion
        seq.reset();
        while (seq.step())
        {
          seq.fillWithIndex(nodeIdx);
          int nn = idx.getIndex(nodeIdx); // index of expansion coefficient
          
// term resulting from orthogonality of Legendre polynomials
          double orthoTerm = 1.0;
          for (unsigned d=0; d<NDIM; ++d)
            orthoTerm *= 2.0/(2*nodeIdx[d]+1.0);
// increment entry
          entry += expandCoeff(nn, k)*expandCoeff(nn, m)*orthoTerm;
        }
        massMatrix(k,m) = entry;
      }
    }
  }

  template <unsigned NDIM>
  void
  LagrangeTensorBasisCalc<NDIM>::calcGradStiff(unsigned dir)
  {
    Lucee::RowMajorIndexer<NDIM> idx(nodeRgn);
    int nodeIdx[NDIM];

// space for matrix
    gradStiff[dir] = Lucee::Matrix<double>(totalNodes, totalNodes);
    gradStiff[dir] = 0.0;

// compute each entry in matrix
    for (unsigned k=0; k<totalNodes; ++k)
    {
      for (unsigned m=0; m<totalNodes; ++m)
      {
        double entry = 0.0;
        for (unsigned r=0; (2*r+1)<nodeRgn.getUpper(dir); ++r)
        {
// this loop over 'r' is basically evaluating \int P'_n(x) P_m(x) dx
// where P_n(x) are Legendre polynomials as a sum of shifted delta
// functions.

// create smaller region
          Lucee::Region<NDIM, int> smallRgn = nodeRgn.resetBounds(dir, 2*r+1, nodeRgn.getUpper(dir));
          Lucee::RowMajorSequencer<NDIM> seq(smallRgn);

// loop over basis function expansion
          while (seq.step())
          {
            seq.fillWithIndex(nodeIdx);
            int nn1 = idx.getIndex(nodeIdx); // index of expansion coefficient into basis k
            nodeIdx[dir] += -(2*r+1); // delta function is shifted by this much
            int nn2 = idx.getIndex(nodeIdx); // index of expansion coefficient into basis m
            
// term resulting from orthogonality of Legendre polynomials
            double orthoTerm = 1.0;
            for (unsigned d=0; d<NDIM; ++d)
            {
              if (d != dir)
                orthoTerm *= 2.0/(2*nodeIdx[d]+1.0);
            }
// increment entry
            entry += 2.0*expandCoeff(nn1, k)*expandCoeff(nn2, m)*orthoTerm;
          }
        }
        gradStiff[dir](k,m) = entry;
      }
    }
  }

  template <unsigned NDIM>
  inline
  Lucee::RowMajorIndexer<NDIM>
  LagrangeTensorBasisCalc<NDIM>::getIndexer() const
  {
    return Lucee::RowMajorIndexer<NDIM>(nodeRgn);
  }

// instantiations
  template class LagrangeTensorBasisCalc<1>;
  template class LagrangeTensorBasisCalc<2>;
  template class LagrangeTensorBasisCalc<3>;
  template class LagrangeTensorBasisCalc<4>;
  template class LagrangeTensorBasisCalc<5>;
}
