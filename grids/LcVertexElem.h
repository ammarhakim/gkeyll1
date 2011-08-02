/**
 * @file	LcVertexElem.h
 *
 * @brief       Vertex element in unstructured grid.
 */

#ifndef LC_VERTEX_ELEMS_H
#define LC_VERTEX_ELEMS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * A class representing a vertex.
 */
  template <typename REAL>
  class VertexElem
  {
    public:
/**
 * Fill with vertex coordinates.
 *
 * @param xv On output vertex coordinates.
 */
      void fillWithCoordinates(REAL xv[3]) const;

    protected:
/**
 * Create vertex given list of vertex coordinates.
 *
 * @param vc Vertex coordinates.
 */
      VertexElem(const std::vector<REAL>& vc);

/**
 * Increment vertex location by one.
 */
      void incr() const
      {
        currLoc += 3; // 3 as (x,y,z) coordinates are stored
      }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
      bool atEnd() const
      {
        if (currLoc<vcoords.size())
          return false;
        return true;
      }

    private:
/** Reference to vertex coordinates */
      const std::vector<REAL>& vcoords;
/** Current location in vcoords array */
      mutable unsigned currLoc;
  };
}

#endif // LC_VERTEX_ELEMS_H
