/**
 * @file	LcUnstructGridCreator.cpp
 *
 * @brief	Unstructured grid creator class.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcUnstructGridCreator.h>

namespace Lucee
{
  template <typename REAL>
  UnstructGridCreator<REAL>::UnstructGridCreator(unsigned ndim)
    : ndim(ndim), currCell(0)
  {
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumVertices(unsigned nv)
  {
    vc.reset(nv);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumCells(unsigned nc)
  {
    c2v.reset(nc);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setVertex(unsigned iv, double xv[3])
  {
// check if there is enough space to add vertex
    if (vc.vcoords.size() < 3*(iv+1))
    {
      Lucee::Except lce("UnstructGridCreator::setVertex: Vertex number ");
      lce << iv << " can not be added: out of space." << std::endl;
      throw lce;
    }

// find correct location and set coordinates
    unsigned loc = 3*iv;
    for (unsigned i=0; i<3; ++i)
      vc.vcoords[loc+i] = xv[i];
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::addTri(unsigned a, unsigned b, unsigned c)
  {
    if (ndim != 2)
      throw Lucee::Except("UnstructGridCreator::addTri: Can add triangle only in 2D grids");
// append connections to connectivity array
    c2v.indices.push_back(a);
    c2v.indices.push_back(b);
    c2v.indices.push_back(c);
// now set offsets correctly
    c2v.offsets[currCell+1] = c2v.offsets[currCell]+3;
    currCell += 1;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::addQuad(unsigned a, unsigned b, unsigned c, unsigned d)
  {
    if (ndim != 2)
      throw Lucee::Except("UnstructGridCreator::addQuad: Can add quad only in 2D grids");
// append connections to connectivity array
    c2v.indices.push_back(a);
    c2v.indices.push_back(b);
    c2v.indices.push_back(c);
    c2v.indices.push_back(d);
// now set offsets correctly
    c2v.offsets[currCell+1] = c2v.offsets[currCell]+4;
    currCell += 1;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::addTet(unsigned a, unsigned b, unsigned c, unsigned d)
  {
    if (ndim != 3)
      throw Lucee::Except("UnstructGridCreator::addTet: Can add tet only in 3D grids");
// append connections to connectivity array
    c2v.indices.push_back(a);
    c2v.indices.push_back(b);
    c2v.indices.push_back(c);
    c2v.indices.push_back(d);
// now set offsets correctly
    c2v.offsets[currCell+1] = c2v.offsets[currCell]+4;
    currCell += 1;
  }

// instantiations
  template class UnstructGridCreator<float>;
  template class UnstructGridCreator<double>;
}
