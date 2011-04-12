/**
 * @file	LcEdgeElem.h
 *
 * @brief       Edge elements in unstructured grids.
 *
 * @version	$Id$
 */

#ifndef LC_EDGE_ELEM_H
#define LC_EDGE_ELEM_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * A class representing a edge.
 */
  template <typename REAL>
  class EdgeElem
  {
    public:
/**
 * Fill with edge center coordinates.
 *
 * @param xv On output coordinates of edge midpoint.
 */
      void fillWithCoordinates(REAL xv[3]) const;

/**
 * Get length of edge.
 *
 * @return lenght of edge.
 */
      REAL getMeasure() const;

/**
 * Fill with normal (unit vector) to edge. Note that this only makes
 * sense when the edge lives in a 2D mesh.
 *
 * @param norm On output normal to edge.
 */
      void fillWithNormal(REAL norm[3]) const;

/**
 * Fill with a tangent (unit) to edge. This particular tangent is in
 * the direction of the edge.
 *
 * @param tng On output tangent to edge.
 */
      void fillWithTangent1(REAL tng[3]) const;

/**
 * Fill with a tangent (unit) to edge. This particular tangent make
 * sense only in 2D and lies in the perpendicular plane. Note that the
 * relation norm = tangent1 X tangent 2 holds.
 *
 * @param tng On output tangent to edge.
 */
      void fillWithTangent2(REAL tng[3]) const;
  };
}

#endif // LC_EDGE_ELEM_H
