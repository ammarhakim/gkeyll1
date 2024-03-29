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
#include <LcColMajorIndexer.h>

// std includes
#include <cmath>

namespace Lucee
{
// set constructor name
  template <> const char *CartProdDecompRegionCalc<1>::id = "CartProd";
  template <> const char *CartProdDecompRegionCalc<2>::id = "CartProd";
  template <> const char *CartProdDecompRegionCalc<3>::id = "CartProd";
  template <> const char *CartProdDecompRegionCalc<4>::id = "CartProd";
  template <> const char *CartProdDecompRegionCalc<5>::id = "CartProd";
  template <> const char *CartProdDecompRegionCalc<6>::id = "CartProd";

  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc()
    : cutIndexer(&Lucee::FixedVector<NDIM, unsigned>( (unsigned) 0)[0], &Lucee::FixedVector<NDIM, int>((int) 1)[0])
  {
// set some reasonable defaults
    unsigned c[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      c[i] = 1;
    setCuts(c);
  }

  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc(const unsigned c[NDIM])
    : cutIndexer(&Lucee::FixedVector<NDIM, unsigned>( (unsigned) 0)[0], &Lucee::FixedVector<NDIM, int>((int) 1)[0])
  {
    setCuts(c);
  }

  template <unsigned NDIM>
  CartProdDecompRegionCalc<NDIM>::CartProdDecompRegionCalc(const int c[NDIM])
    : cutIndexer(&Lucee::FixedVector<NDIM, unsigned>( (unsigned) 0)[0], &Lucee::FixedVector<NDIM, int>((int) 1)[0])
  {
    unsigned cc[NDIM];
    for (unsigned d=0; d<NDIM; ++d) cc[d] = c[d];
    setCuts(cc);
  }

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    DecompRegionCalcIfc<NDIM>::readInput(tbl);

    if (tbl.hasNumVec("cuts"))
    {
      std::vector<double> dblCuts = tbl.getNumVec("cuts");
      if (dblCuts.size() != NDIM)
      {
        Lucee::Except lce("CartProdDecompRegionCalc::readInput: 'cuts' table must have ");
        lce << NDIM << " entries. Only " << dblCuts.size() << " provided";
        throw lce;
      }
      unsigned myCuts[NDIM];
      for (unsigned i=0; i<NDIM; ++i)
        myCuts[i] = (unsigned) dblCuts[i];
      setCuts(myCuts);
    }
    else
    {
      throw Lucee::Except("CartProdDecompRegionCalc::readInput: Must provide 'cuts' table.");
    }
  }

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::fillWithCuts(int c[NDIM]) const
  {
    for (unsigned d=0; d<NDIM; ++d)
      c[d] = cuts[d];
  }

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn)
  {
    if (nrgs != nsub)
    {
      Lucee::Except lce("CartProdDecompRegionCalc::decompose: Number of sub-regions from cuts (");
      lce << nsub << ") does not match requested sub-regions (" << nrgs << ")";
      throw lce;
    }

// struct to hold shapes along each direction
    struct { std::vector<unsigned> shape, start; } dimShapes[NDIM];

// NOTE: The basic idea is to first set shape[dim] to be the nearest
// smaller integer number of cells. Then, the remainder cells from the
// division are added one by one to each of the cuts. This ensures
// that each cut differs only by a single cell in each dim.

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

  template <unsigned NDIM>
  void
  CartProdDecompRegionCalc<NDIM>::setCuts(const unsigned c[NDIM])
  {
    nsub = 1;
    int zeros[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      cuts[i] = c[i];
      nsub *= cuts[i];
      zeros[i] = 0;
    }
// reset indexer
    cutIndexer = Lucee::ColMajorIndexer<NDIM>(c, zeros);
  }

// instantiations
  template class CartProdDecompRegionCalc<1>;
  template class CartProdDecompRegionCalc<2>;
  template class CartProdDecompRegionCalc<3>;
  template class CartProdDecompRegionCalc<4>;
  template class CartProdDecompRegionCalc<5>;
  template class CartProdDecompRegionCalc<6>;
}
