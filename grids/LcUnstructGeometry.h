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
// forward declare creator class
  template <typename REAL> class UnstructGridCreator;
// forward declare grid class
  template <typename REAL> class UnstructGrid;

/**
 * Class to hold unstructured grid geometry information. This class is
 * private and hence can not be accessed directly.
 */
  template <unsigned NDIM, typename REAL>
  class UnstructGeometry
  {
// declare friends
      template <typename RT> friend class UnstructGridCreator;
      template <typename RT> friend class UnstructGrid;

/**
 * Create empty geometery object. The 'reset' method must be called to
 * allocate any memory to store the geometry.
 */
      UnstructGeometry();

/** 
 * Reset geometry object to contain specified number of vertices.
 *
 * @param nv Number of vertices.
 */
      void reset(unsigned nv);

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
