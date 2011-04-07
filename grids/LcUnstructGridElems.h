/**
 * @file	LcUnstructGridElems.h
 *
 * @brief File to define basic elements in unstructured meshes.
 *
 * @version	$Id$
 */

#ifndef LC_UNSTRUCT_GRID_ELEMS_H
#define LC_UNSTRUCT_GRID_ELEMS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
// declare grid class
  template <typename REAL> class UnstructGrid;

/**
 * A class representig a vertex.
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
 * @param vc List of vertex coordinates.
 */
      VertexElem(const std::vector<REAL>& vc);

/**
 * Increment vertex location by one.
 */
      void incr() const;

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
      std::vector<REAL> vcoords;
/** Current location in vcoords array */
      mutable unsigned currLoc;
  };

/**
 * A basic template class that allows accessing the grid elements
 * using template parameter.
 */
  template <typename REAL, unsigned NDIM> class GridElem;

/**
 * Wrapper for 0d elements (vertex).
 */
  template <typename REAL>
  class GridElem<REAL, 0> : public VertexElem<REAL>
  {
// declare grid friend so it can fiddle around with privates
      template <typename R> friend class UnstructGrid;
    private:
/**
 * Initialize element.
 *
 * @param vc List of vertex coordinates.
 */
      GridElem(const std::vector<REAL>& vc)
        : VertexElem<REAL>(vc)
      {
      }

/**
 * Increment vertex pointer by one element.
 */
      void incr() const { VertexElem<REAL>::incr(); }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
      bool atEnd() const { return VertexElem<REAL>::atEnd(); }
  };
}

#endif // LC_UNSTRUCT_GRID_ELEMS_H
