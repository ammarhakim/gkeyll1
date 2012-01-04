/**
 * @file	LcCartProdDecompRegionCalc.cpp
 *
 * @brief	Decomposition with specified cuts.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcMatrix.h>

namespace Lucee
{
  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc(const unsigned c[NDIM])
  {
    nsub = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      cuts[i] = c[i];
      nsub *= cuts[i];
    }
  }

  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc(const int c[NDIM])
  {
    nsub = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      cuts[i] = c[i];
      nsub *= cuts[i];
    }
  }

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn)
  {
// check if cuts match up with number of regions
    if (nrgs != nsub)
    {
      Lucee::Except lce("CartProdDecompRegionCalc::decompose: Number of sub-regions from cuts (");
      lce << nsub << ") does not match requested sub-regions (" << nrgs << ")";
      throw lce;
    }

// struct to hold shapes along each direction
    struct { std::vector<unsigned> shape, start; } dimShapes[NDIM];

// compute shapes in each direction
    for (unsigned dim=0; dim<NDIM; ++dim)
    {
      dimShapes[dim].shape.resize(cuts[dim]);
      dimShapes[dim].start.resize(cuts[dim]);

 // number of cells in each sub-region by default
      int baseShape = globalRgn.getShape(dim)/cuts[dim];
// set default shape of each sub-region
      for (unsigned r=0; r<cuts[dim]; ++r)
        dimShapes[dim].shape[r] = baseShape;

// remainder cells
      int extraCells = globalRgn.getShape(dim) % cuts[dim];
      if (extraCells > 0)
      {
// redistribute extra cells to each sub-region
        unsigned r = 0; // first region to add extra cell
        for (unsigned e=0; e<extraCells; ++e)
        {
          dimShapes[dim].shape[r % cuts[dim]] += 1;
          r++;
        }
      }

      dimShapes[dim].start[0] = globalRgn.getLower(dim);
// set start index of each sub-region
      for (unsigned r=1; r<cuts[dim]; ++r)
        dimShapes[dim].start[r] = 
          dimShapes[dim].start[r-1] + dimShapes[dim].shape[r-1];
    }

    int idx[NDIM];
    Lucee::Region<NDIM, int> cutsRgn(cuts);
// add regions to decomposition
    typename Lucee::ColMajorSequencer<NDIM> seq(cutsRgn);
    while (seq.step())
    { // loop over region defined by cuts
      seq.fillWithIndex(idx);
      int lower[NDIM], upper[NDIM];
      for (unsigned d=0; d<NDIM; ++d)
      {
        lower[d] = dimShapes[d].start[idx[d]];
        upper[d] = lower[d] + dimShapes[d].shape[idx[d]];
      }
      this->addRegion(Lucee::Region<NDIM, int>(lower, upper));
    }
  }

// instantiations
  template class CartProdDecompRegionCalc<1>;
  template class CartProdDecompRegionCalc<2>;
  template class CartProdDecompRegionCalc<3>;
  template class CartProdDecompRegionCalc<4>;
  template class CartProdDecompRegionCalc<5>;
  template class CartProdDecompRegionCalc<6>;
  template class CartProdDecompRegionCalc<7>;
}
