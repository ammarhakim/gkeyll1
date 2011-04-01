/**
 * @file	LcUnstructGeometry.h
 *
 * @brief	Class holding geometry of unstructured grids.
 *
 * @version	$Id$
 */

#ifndef LC_UNSTRUCT_GEOMETRY_H
#define LC_UNSTRUCT_GEOMETRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to hold unstructured grid geometry information.
 */
  template <unsigned NDIM, typename REAL>
  class UnstructGeometry
  {
/** 
 * Create a geometry object with specified number of vertices.
 *
 * @param nv Number of vertices.
 */
      UnstructGeometry(unsigned nv);

    private:
/** Vertex coordinates NDIM*numVertices stored as x,y,z */
      std::vector<REAL> vcoords;

/** No copying allowed */
      UnstructGeometry(const UnstructGeometry<NDIM, REAL>&);
/** No assignment allowed */
      UnstructGeometry<NDIM, REAL>& operator=(const UnstructGeometry<NDIM, REAL>&);
  };
}

#endif //  LC_UNSTRUCT_GEOMETRY_H
