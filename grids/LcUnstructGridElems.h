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

// lucee includes
#include <LcCellElem.h>
#include <LcEdgeElem.h>
#include <LcUnstructGeometry.h>
#include <LcVertexElem.h>

// std includes
#include <vector>

namespace Lucee
{
// declare grid class
  template <typename REAL> class UnstructGrid;

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
 * @param geom Geometry object.
 */
      GridElem(const Lucee::UnstructGeometry<3, REAL>& geom)
        : VertexElem<REAL>(geom.vcoords)
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
/**
 * Initialize element.
 *
 * @param geom Geometry object.
 */
      GridElem(const Lucee::UnstructGeometry<3, REAL>& geom)
        : CellElem<REAL>(geom.cellCentroid, geom.cellVolume)
      {
      }

/**
 * Increment vertex pointer by one element.
 */
      void incr() const { CellElem<REAL>::incr(); }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
      bool atEnd() const { return CellElem<REAL>::atEnd(); }
  };
}

#endif // LC_UNSTRUCT_GRID_ELEMS_H
