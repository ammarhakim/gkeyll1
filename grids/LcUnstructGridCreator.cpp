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
  static short tet_c = TET_CELL_T;
  static short hex_c = HEX_CELL_T;
  static short tri_c = TRI_CELL_T;
  static short quad_c = QUAD_CELL_T;

  template <typename REAL>
  UnstructGridCreator<REAL>::UnstructGridCreator(unsigned ndim)
    : ndim(ndim), currCell(0)
  {
// initialize cell count
    cellCount[tet_c] = 0;
    cellCount[hex_c] = 0;
    cellCount[tri_c] = 0;
    cellCount[quad_c] = 0;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::fillWithGeometry(Lucee::UnstructGeometry<3, REAL>& geo) const
  {
    geo = vc;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::fillWithConnectivity(Lucee::UnstructConnectivity& conn) const
  {
    conn = c2v;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::fillWithCellType(std::vector<short>& ctype) const
  {
    ctype = cellType;
  }

  template <typename REAL>
  unsigned
  UnstructGridCreator<REAL>::getNumTri() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TRI_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGridCreator<REAL>::getNumQuad() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(QUAD_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGridCreator<REAL>::getNumTet() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(TET_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  unsigned
  UnstructGridCreator<REAL>::getNumHex() const
  {
    std::map<short, unsigned>::const_iterator itr
      = cellCount.find(HEX_CELL_T);
    return itr->second;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumVertices(unsigned nv)
  {
    vc.setNumVertices(nv);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setNumCells(unsigned nc)
  {
    c2v.reset(nc);
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setVertex(unsigned iv, REAL xv[3])
  {
// check if there is enough space to add vertex
    if (vc.vcoords.size() < 3*(iv+1))
    {
      Lucee::Except lce("UnstructGridCreator::setVertex: Vertex number ");
      lce << iv << " can not be added: out of space." << std::endl;
      throw lce;
    }
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setVertexXCoord(unsigned iv, REAL x)
  {
// check if there is enough space to add vertex
    if (vc.vcoords.size() < 3*(iv+1))
    {
      Lucee::Except lce("UnstructGridCreator::setVertex: Vertex number ");
      lce << iv << " can not be added: out of space." << std::endl;
      throw lce;
    }

// find correct location and set coordinates
    unsigned loc = 3*iv+0;
    vc.vcoords[loc] = x;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setVertexYCoord(unsigned iv, REAL x)
  {
// check if there is enough space to add vertex
    if (vc.vcoords.size() < 3*(iv+1))
    {
      Lucee::Except lce("UnstructGridCreator::setVertex: Vertex number ");
      lce << iv << " can not be added: out of space." << std::endl;
      throw lce;
    }

// find correct location and set coordinates
    unsigned loc = 3*iv+1;
    vc.vcoords[loc] = x;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::setVertexZCoord(unsigned iv, REAL x)
  {
// check if there is enough space to add vertex
    if (vc.vcoords.size() < 3*(iv+1))
    {
      Lucee::Except lce("UnstructGridCreator::setVertex: Vertex number ");
      lce << iv << " can not be added: out of space." << std::endl;
      throw lce;
    }

// find correct location and set coordinates
    unsigned loc = 3*iv+2;
    vc.vcoords[loc] = x;
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
// set cell type
    cellType.push_back(tri_c);
    currCell += 1;
    cellCount[tri_c] += 1;
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
    cellType.push_back(quad_c);
    currCell += 1;
    cellCount[quad_c] += 1;
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
    cellType.push_back(tet_c);
    currCell += 1;
    cellCount[tet_c] += 1;
  }

  template <typename REAL>
  void
  UnstructGridCreator<REAL>::addHex(unsigned a1, unsigned b1, unsigned c1, unsigned d1,
    unsigned a2, unsigned b2, unsigned c2, unsigned d2)
  {
    if (ndim != 3)
      throw Lucee::Except("UnstructGridCreator::addHex: Can add hex only in 3D grids");
// append connections to connectivity array
    c2v.indices.push_back(a1);
    c2v.indices.push_back(b1);
    c2v.indices.push_back(c1);
    c2v.indices.push_back(d1);
    c2v.indices.push_back(a2);
    c2v.indices.push_back(b2);
    c2v.indices.push_back(c2);
    c2v.indices.push_back(d2);
// now set offsets correctly
    c2v.offsets[currCell+1] = c2v.offsets[currCell]+8;
    cellType.push_back(hex_c);
    currCell += 1;
    cellCount[hex_c] += 1;
  }

// instantiations
  template class UnstructGridCreator<float>;
  template class UnstructGridCreator<double>;
}
