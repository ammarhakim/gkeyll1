/**
 * @file	LcDecompRegionCalcIfc.cpp
 *
 * @brief	Base class for decomposition algorithms.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>
#include <LcRowMajorSequencer.h>

// std includes
#include <iostream>

namespace Lucee
{
// set module name
  template <> const char *DecompRegionCalcIfc<1>::id = "DecompRegionCalc1D";
  template <> const char *DecompRegionCalcIfc<2>::id = "DecompRegionCalc2D";
  template <> const char *DecompRegionCalcIfc<3>::id = "DecompRegionCalc3D";
  template <> const char *DecompRegionCalcIfc<4>::id = "DecompRegionCalc4D";
  template <> const char *DecompRegionCalcIfc<5>::id = "DecompRegionCalc5D";

  template <unsigned NDIM> 
  DecompRegionCalcIfc<NDIM>::DecompRegionCalcIfc()
  {
    for (unsigned d=0; d<NDIM; ++d)
      isPeriodic[d] = false;
  }

  template <unsigned NDIM> 
  void
  DecompRegionCalcIfc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    for (unsigned d=0; d<NDIM; ++d)
      isPeriodic[d] = false;
// check if any directions are periodic
    if (tbl.hasNumVec("periodicDirs"))
    {
      std::vector<double> pd = tbl.getNumVec("periodicDirs");
      for (unsigned i=0; i<pd.size(); ++i)
        setPeriodicDir( (unsigned) pd[i]);
    }
  }

  template <unsigned NDIM> 
  void
  DecompRegionCalcIfc<NDIM>::setPeriodicDir(unsigned dir)
  {
    if (dir>=0 && dir<NDIM)
      isPeriodic[dir] = true;
    else
    {
      Lucee::Except lce(
        "DecompRegionCalcIfc::setPeriodicDirs: Direction must be between 0 and ");
      lce << (NDIM-1);
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::calcDecomp(unsigned nrgns, Lucee::DecompRegion<NDIM>& decompRgn)
  {
    decompRgnPtr = &decompRgn; // so addRegion can access
// clear existing decomposition
    decompRgn.clearDecomp();
// call derived class to compute decomposition
    this->decompose(nrgns, decompRgn.getGlobalRegion());
// check if decompostion covers global region
    if (decompRgn.checkCovering() == false)
      throw Lucee::Except(
        "DecompRegionCalcIfc: Decomposition does not cover global region");
// take into account periodic directions (if any)
    handlePeriodicDirs();
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::addRegion(const Lucee::Region<NDIM, int>& subRgn)
  {
    decompRgnPtr->addRegion(subRgn);
  }

  template <unsigned NDIM>
  void
  DecompRegionCalcIfc<NDIM>::handlePeriodicDirs()
  {
    Lucee::Region<NDIM, int> rgn = decompRgnPtr->getGlobalRegion();
// make a copy
    std::vector<Lucee::Region<NDIM, int> > realVec;
    for (unsigned i=0; i<decompRgnPtr->getNumRegions(); ++i)
      realVec.push_back( decompRgnPtr->getRegion(i) );

    // for (typename BoxMap_t::const_iterator i = realMap.begin(); i != realMap.end(); ++i) 
    // {
    //   BoxRPair_t insertme((*i).first, (*i).first);
    //   data->boxRank.insert(insertme);
    // }

    int itrMin[NDIM], itrMax[NDIM];
    for (size_t i = 0; i < NDIM; ++i) {
      itrMin[i] = isPeriodic[i] ? -1 : 0;
      itrMax[i] = isPeriodic[i] ? 2 : 1;
    }
    Lucee::Region<NDIM, int> shift(itrMin, itrMax);
    Lucee::RowMajorSequencer<NDIM> seq(shift);
    size_t nTotBoxes = decompRgnPtr->getNumRegions();
    int idx[NDIM];
    while ( seq.step() ) {
      seq.fillWithIndex(idx);
// skip case of zero offset
      bool skipStep = true;
      for (size_t i = 0; i < NDIM; ++i) {
        if (idx[i] != 0) {
          skipStep = false;
          break;
        }
      }
// shift boxes, renumber, and add to the box map
      if (!skipStep) {
        int lext[NDIM], uext[NDIM], dist[NDIM];
        for (size_t i = 0; i < NDIM; ++i) {
          dist[i] = rgn.getShape(i);
          lext[i] = -1*idx[i]*dist[i];
          uext[i] =  1*idx[i]*dist[i];
        }
        typename std::vector< Lucee::Region<NDIM, int> >::const_iterator itr
          = realVec.begin();
        for ( ; itr != realVec.end(); ++itr) {
          //std::cout << "nTotBoxes " << nTotBoxes << std::endl;
          Lucee::Region<NDIM, int> rgn = itr->extend(lext, uext);

          // for (unsigned dd=0; dd<NDIM; ++dd)
          //   std::cout << rgn.getLower(dd) << ", ";
          // std::cout << std::endl;
          // for (unsigned dd=0; dd<NDIM; ++dd)
          //   std::cout << rgn.getUpper(dd) << ", ";
          // std::cout << std::endl;
          
          // BoxPair_t insertme(nTotBoxes, (*i).second.extend(lext, uext));
          // BoxRPair_t insertmetoo(nTotBoxes, (*i).first);
          // data->boxMap.insert(insertme);
          // data->boxRank.insert(insertmetoo);
          nTotBoxes++;
        }
      }
    }
  }

// instantiations
  template class DecompRegionCalcIfc<1>;
  template class DecompRegionCalcIfc<2>;
  template class DecompRegionCalcIfc<3>;
  template class DecompRegionCalcIfc<4>;
  template class DecompRegionCalcIfc<5>;
}
