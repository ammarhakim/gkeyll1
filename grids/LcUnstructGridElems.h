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
 * @param vc List of vertex coordinates.
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
      std::vector<REAL> vcoords;
/** Current location in vcoords array */
      mutable unsigned currLoc;
  };

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

/**
 * A class representing a face.
 */
  template <typename REAL>
  class FaceElem
  {
    public:
/**
 * Fill with face centroid coordinates.
 *
 * @param xv On output coordinates of face centroid.
 */
      void fillWithCoordinates(REAL xv[3]) const;

/**
 * Get area of face.
 *
 * @return area of face.
 */
      REAL getMeasure() const;

/**
 * Fill with normal (unit vector) to face.
 *
 * @param norm On output normal to face.
 */
      void fillWithNormal(REAL norm[3]) const;

/**
 * Fill with a tangent (unit) to face.
 *
 * @param tng On output tangent to face.
 */
      void fillWithTangent1(REAL tng[3]) const;

/**
 * Fill with a tangent (unit) to face. Note that the relation norm =
 * tangent1 X tangent 2 holds.
 *
 * @param tng On output tangent to face.
 */
      void fillWithTangent2(REAL tng[3]) const;
  };

/**
 * A class representing a cell.
 */
  template <typename REAL>
  class CellElem
  {
    public:
/**
 * Fill with cell centroid coordinates.
 *
 * @param xv On output coordinates of cell centroid.
 */
      void fillWithCoordinates(REAL xv[3]) const;

/**
 * Get volume of cell.
 *
 * @return volume of cell.
 */
      REAL getMeasure() const;
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

/**
 * Wrapper for 1d elements (edge).
 */
  template <typename REAL>
  class GridElem<REAL, 1> : public EdgeElem<REAL>
  {
// declare grid friend so it can fiddle around with privates
      template <typename R> friend class UnstructGrid;
    private:
  };

/**
 * Wrapper for 2d elements (face).
 */
  template <typename REAL>
  class GridElem<REAL, 2> : public FaceElem<REAL>
  {
// declare grid friend so it can fiddle around with privates
      template <typename R> friend class UnstructGrid;
    private:
  };

/**
 * Wrapper for 3d elements (cell).
 */
  template <typename REAL>
  class GridElem<REAL, 3> : public CellElem<REAL>
  {
// declare grid friend so it can fiddle around with privates
      template <typename R> friend class UnstructGrid;
    private:
  };
}

#endif // LC_UNSTRUCT_GRID_ELEMS_H
