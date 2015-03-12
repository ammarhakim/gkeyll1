/**
 * @file	LcSubCartProdDecompRegionCalc.cpp
 *
 * @brief	Decomposition based on higher-dimensional decomposition
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcSubCartProdDecompRegionCalc.h>

namespace Lucee
{
// set constructor names: this seems rather strange, but essentially,
// we are adding a set of constructors to each of the 1D, 2D, 3D
// decomposition calculator modules. We stop at 3D as this is the
// maximum dimension of the configuration space we expect in kinetic
// simulations.

  template <> const char *SubCartProdDecompRegionCalc<1,2>::id = "SubCartProd2D";
  template <> const char *SubCartProdDecompRegionCalc<1,3>::id = "SubCartProd3D";
  template <> const char *SubCartProdDecompRegionCalc<1,4>::id = "SubCartProd4D";
  template <> const char *SubCartProdDecompRegionCalc<1,5>::id = "SubCartProd5D";
  template <> const char *SubCartProdDecompRegionCalc<1,6>::id = "SubCartProd6D";

  template <> const char *SubCartProdDecompRegionCalc<2,3>::id = "SubCartProd3D";
  template <> const char *SubCartProdDecompRegionCalc<2,4>::id = "SubCartProd4D";
  template <> const char *SubCartProdDecompRegionCalc<2,5>::id = "SubCartProd5D";
  template <> const char *SubCartProdDecompRegionCalc<2,6>::id = "SubCartProd6D";

  template <> const char *SubCartProdDecompRegionCalc<3,4>::id = "SubCartProd4D";
  template <> const char *SubCartProdDecompRegionCalc<3,5>::id = "SubCartProd5D";
  template <> const char *SubCartProdDecompRegionCalc<3,6>::id = "SubCartProd6D";

/**
 * Check if communicator is valid. This basically works under the
 * assumption that an invalid communicator is NULL.
 *
 * @param c Communicator to check.
 * @return true if communicator is valid, false otherwise.
 */
  static bool _isValid(TxCommBase *c) { return (bool) c; }

  template <unsigned NDIM, unsigned HDIM>
  SubCartProdDecompRegionCalc<NDIM, HDIM>::SubCartProdDecompRegionCalc()
  {
// set some reasonable defaults
    int c[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      c[i] = 1;
    setCuts(c);
  }

  template <unsigned NDIM, unsigned HDIM>
  void
  SubCartProdDecompRegionCalc<NDIM, HDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    DecompRegionCalcIfc<NDIM>::readInput(tbl);
// get reference to parent decomposition
    typename Lucee::CartProdDecompRegionCalc<HDIM>& decompCalc
        = tbl.template getObject<Lucee::CartProdDecompRegionCalc<HDIM> >("decomposition");

// parent cuts
    int pCuts[HDIM];
    decompCalc.fillWithCuts(pCuts);
// parent communicator
    TxCommBase *pComm = decompCalc.getComm();
    
// list of directions to collect
    std::vector<double> fdc = tbl.getNumVec("collectDirections");
    if (fdc.size() != NDIM)
    {
      Lucee::Except lce("SubCartProdDecompRegionCalc::readInput: 'collectDirections' should have exactly ");
      lce << NDIM << " entries. Specified " << fdc.size() << " instead";
      throw lce;
    }
    std::vector<unsigned> dc(NDIM);
    for (unsigned d=0; d<NDIM; ++d) dc[d] = (unsigned) fdc[d];

// determine cuts for sub-decomposition
    int cCuts[NDIM];
    for (unsigned d=0; d<NDIM; ++d) cCuts[d] = pCuts[dc[d]];
    setCuts(cCuts);

    Lucee::Region<HDIM, int> pCutRgn(pCuts); // parent region
    Lucee::Region<NDIM, int> cCutRgn(cCuts); // child-decomp region

    Lucee::Region<HDIM, int> defRgn(pCutRgn);
// create new region with collected dimensions deflated
    for (unsigned d=0; d<NDIM; ++d)
      defRgn = defRgn.deflate(dc[d]);

    unsigned numSubComms = defRgn.getVolume();

// create two sequencers, one over deflated region, and other over
// collected directions
    Lucee::ColMajorSequencer<HDIM> defSeq(defRgn);
    Lucee::ColMajorSequencer<NDIM> colSeq(cCutRgn);
// indexer over full parent cuts region
    Lucee::ColMajorIndexer<HDIM> pIdxr(pCutRgn);

// now build list of child communicators
    std::vector<TxCommBase*> subComms;

    int defIdx[HDIM], colIdx[NDIM];
    unsigned count = 0;
    while (defSeq.step())
    {
      defSeq.fillWithIndex(defIdx);

      std::vector<int> subRanks;
      colSeq.reset();
      while (colSeq.step())
      {
        colSeq.fillWithIndex(colIdx);
        for (unsigned d=0; d<NDIM; ++d) defIdx[dc[d]] = colIdx[d];
// This works because CartProdDecompRegionCalc class uses column major
// order to index region decomposition. If that were to change, all
// hell would break loose. It may be best to rationalize use of a
// particular ordering throughout. (Ammar Hakim: 3/11/2015).
        subRanks.push_back(pIdxr.getIndex(defIdx));
      }
// comm-split parent communicator
      subComms.push_back(pComm->createSubComm(subRanks));
    }

// finally: set valid communicator for this decomposition
    setValidComm(subComms);
  }

  template <unsigned NDIM, unsigned HDIM>
  void
  SubCartProdDecompRegionCalc<NDIM, HDIM>::decompose(unsigned nrgs, const Lucee::Region<NDIM, int>& globalRgn)
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

  template <unsigned NDIM, unsigned HDIM>
  void
  SubCartProdDecompRegionCalc<NDIM, HDIM>::setCuts(const int c[NDIM])
  {
    nsub = 1;
    for (unsigned i=0; i<NDIM; ++i)
    {
      cuts[i] = c[i];
      nsub *= cuts[i];
    }
  }

  template <unsigned NDIM, unsigned HDIM>
  void
  SubCartProdDecompRegionCalc<NDIM, HDIM>::setValidComm(const std::vector<TxCommBase*>& clist)
  {
    for (unsigned i=0; i<clist.size(); ++i)
      if (_isValid(clist[i]))
        this->setComm(clist[i]);
  }

// instantiations
  template class SubCartProdDecompRegionCalc<1,2>;
  template class SubCartProdDecompRegionCalc<1,3>;
  template class SubCartProdDecompRegionCalc<1,4>;
  template class SubCartProdDecompRegionCalc<1,5>;
  template class SubCartProdDecompRegionCalc<1,6>;

  template class SubCartProdDecompRegionCalc<2,3>;
  template class SubCartProdDecompRegionCalc<2,4>;
  template class SubCartProdDecompRegionCalc<2,5>;
  template class SubCartProdDecompRegionCalc<2,6>;

  template class SubCartProdDecompRegionCalc<3,4>;
  template class SubCartProdDecompRegionCalc<3,5>;
  template class SubCartProdDecompRegionCalc<3,6>;
}
