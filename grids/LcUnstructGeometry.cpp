/**
 * @file	LcUnstructGeometry.cpp
 *
 * @brief	Class holding geometry of unstructured grids.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructGeometry.h>

namespace Lucee
{
  template <unsigned NDIM, typename REAL>
  UnstructGeometry<NDIM, REAL>::UnstructGeometry()
  {
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumVertices(unsigned nv)
  {
    vcoords.clear();
    vcoords.resize(NDIM*nv);
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumCells(unsigned nc)
  {
    cellCentroid.clear();
    cellCentroid.resize(NDIM*nc);
    cellVolume.clear();
    cellVolume.resize(nc);
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumFaces(unsigned nf, bool storeNormal, bool storeTangents)
  {
    faceCenter.clear();
    faceCenter.resize(NDIM*nf);
    faceArea.clear();
    faceArea.resize(nf);
// allocate normal and tangents if requested
    if (storeNormal) 
    {
      faceNormal.clear();
      faceNormal.resize(NDIM*nf);
    }
    if (storeTangents) 
    {
      faceTangent1.clear();
      faceTangent1.resize(NDIM*nf);
      faceTangent2.clear();
      faceTangent2.resize(NDIM*nf);
    }
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumEdges(unsigned nf, bool storeNormal, bool storeTangents)
  {
    edgeCenter.clear();
    edgeCenter.resize(NDIM*nf);
    edgeLength.clear();
    edgeLength.resize(nf);
// allocate normal and tangents if requested
    if (storeNormal) 
    {
      edgeNormal.clear();
      edgeNormal.resize(NDIM*nf);
    }
    if (storeTangents) 
    {
      edgeTangent1.clear();
      edgeTangent1.resize(NDIM*nf);
      edgeTangent2.clear();
      edgeTangent2.resize(NDIM*nf);
    }
  }

// instantiations
  template class Lucee::UnstructGeometry<1, float>;
  template class Lucee::UnstructGeometry<2, float>;
  template class Lucee::UnstructGeometry<3, float>;

  template class Lucee::UnstructGeometry<1, double>;
  template class Lucee::UnstructGeometry<2, double>;
  template class Lucee::UnstructGeometry<3, double>;
}
