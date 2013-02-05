/**
 * @file	LcCartGeneralDecompRegionCalc.cpp
 *
 * @brief	Decomposition with specified cuts.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCartGeneralDecompRegionCalc.h>
#include <LcMatrix.h>

// std includes
#include <cmath>

namespace Lucee
{
// set constructor name
  template <> const char *CartGeneralDecompRegionCalc<1>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<2>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<3>::id = "CartGeneral";
  template <> const char *CartGeneralDecompRegionCalc<4>::id = "CartGeneral";

  template <unsigned NDIM>
  CartGeneralDecompRegionCalc<NDIM>::CartGeneralDecompRegionCalc()
  {
  }

  template <unsigned NDIM>
  void
  CartGeneralDecompRegionCalc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    DecompRegionCalcIfc<NDIM>::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  CartGeneralDecompRegionCalc<NDIM>::decompose(unsigned num, const Lucee::Region<NDIM, int>& decompBox)
  {
// This algorithm was originally written by Mahmood Miah and A. Hakim
// for the Facets project. This implementation is essentially the same
// as in Facets. I am retaining the original comments. (A. Hakim,
// 2/5/2013)

//  begin by determining the number of breaks in the first N-1 dimension
    double side = std::pow((double) num, (double) 1 / NDIM);
    size_t sbounds[2];
    sbounds[0] = (size_t) std::floor(side);
    sbounds[1] = (size_t) std::ceil(side);
    size_t boundtest = 1;
    size_t breaks[NDIM];
    for (size_t i = 0; i < NDIM; ++i) {
      boundtest *= sbounds[0];
    }
//  This does a series of comparisons to figure out how to distribute
//  boxes in the first N-1 dimensions
    for (size_t i = 0; i < NDIM; ++i) {
      boundtest = boundtest * sbounds[1] / sbounds[0];
      if (num <= boundtest) {
        for (size_t j = 0; j < i; ++j) {
          breaks[j] = sbounds[1];
        }
        for (size_t j = i; (int) j < (int) NDIM - 1; ++j) {
          breaks[j] = sbounds[0];
        }
        break;
      }
    }
    breaks[NDIM-1] = 1;

//  Now break down the box into a series of strips
    std::vector< Lucee::Region<NDIM, int> > intermediate = breakBoxes(breaks, decompBox);

// Do the final breakdown, the first rmdr boxes have sbrk + 1 breaks
// in it, while the rest have sbrk breaks in it.  rmdr is the
// remainder of num strips into num processes.  sbrk is the quotient
// of num strips into num processes.  Also rebalance boxes.
    size_t rmdr = num % intermediate.size();
    size_t sbrk = num / intermediate.size();
    size_t dims1[NDIM];
    std::vector<size_t> dimsary;
    for (size_t i = 0; (int) i < (int) NDIM - 1; ++i) {
      dims1[i] = breaks[i];
    }
//  This whole godawful messy thing rebalances boxes to all have ~
//  equal areas.
    dimsary.resize(intermediate.size());
    for (size_t i = 0; i < intermediate.size(); ++i) {
      dimsary[i] = i < rmdr ? sbrk + 1 : sbrk;
    }
    for (size_t i = 0; (int) i < (int) NDIM - 1; ++i) {
      size_t step = 1;
      for (size_t j = 0; j < NDIM - 1 - i; ++j) {
        step *= dims1[j];
      }
      size_t substep = 1;
      for (int j = 0; j < (int)NDIM - 2 - (int)i; ++j) {
        substep *= dims1[j];
      }
      for (size_t n = 0; n < intermediate.size(); n += step) {
        size_t weightsum = 0;
        for (size_t j = n; j < n + step; ++j) {
          weightsum += dimsary[j];
        }
        size_t rollingsum = 0;
        for (size_t j = n; j < n + step; j += substep) {
          size_t stepsum = 0;
          for (size_t k = j; k < j + substep; ++k) {
            stepsum += dimsary[k];
          }
          rollingsum += stepsum;
          int newUpper = rollingsum * decompBox.getShape(NDIM-2-i) / weightsum;
          for (size_t k = j; k < j + substep; ++k) {
            intermediate[k].setUpper(NDIM-2-i, newUpper);
            if ( j + substep < n + step ) {
              intermediate[k+substep].setLower(NDIM-2-i, newUpper);
            }
          }
        }
      }
    }
    for (size_t i = 0; (int) i < (int) NDIM - 1; ++i) {
      breaks[i] = 1;
    }
    size_t count = 0;
    for (size_t i = 0; i < intermediate.size(); ++i) {
      breaks[NDIM-1] = i < rmdr ? sbrk + 1 : sbrk;
      std::vector< Lucee::Region<NDIM, int> > temp;
      temp = breakBoxes(breaks, intermediate[i]);
// insert boxes
      typename std::vector<Lucee::Region<NDIM, int> >::const_iterator j;
      for (j = temp.begin(); j != temp.end(); ++j)
        this->addRegion( *j );
    }
  }

  template <unsigned NDIM>
  std::vector<Lucee::Region<NDIM, int> >
  CartGeneralDecompRegionCalc<NDIM>::breakBoxes(size_t subBox[NDIM], 
    const Lucee::Region<NDIM, int>& box) {
    std::vector<Lucee::Region<NDIM, int> > boxes;
 // compute total number of boxes
    size_t nbrks = 1;
    for (size_t i=0; i<NDIM; ++i)
      nbrks *= subBox[i];
// construct boxes
    size_t itr[NDIM];
    for (size_t i=0; i<NDIM; ++i)
      itr[i] = 0;
    int newLower[NDIM], newUpper[NDIM];
// COULD REPLACE THIS WITH A SEQUENCER OBJECT
    for (size_t i=0; i<nbrks; ++i) {
      for (size_t j=0; j<NDIM; ++j) {
        newLower[j] = box.getLower(j) + itr[j] * box.getShape(j) / subBox[j];
        newUpper[j] = box.getLower(j) + (itr[j] + 1) * box.getShape(j) / subBox[j];
      }
      boxes.push_back(Lucee::Region<NDIM, int>(newLower, newUpper));

      for (size_t j=0; j<NDIM; ++j) {
        if (++itr[j] >= subBox[j])
          itr[j] = 0;
        else
          break;
      }
    }
    return boxes;
  }

// instantiations
  template class CartGeneralDecompRegionCalc<1>;
  template class CartGeneralDecompRegionCalc<2>;
  template class CartGeneralDecompRegionCalc<3>;
  template class CartGeneralDecompRegionCalc<4>;
}
