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
#include <LcFaceElem.h>
#include <LcGridGeometry.h>
#include <LcVertexElem.h>

// std includes
#include <vector>

namespace Lucee
{
// declare grid class
  template <typename REAL> class UnstructGrid;

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
      GridElem(const Lucee::GridGeometry<3, REAL>& geom)
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
/**
 * Initialize element.
 *
 * @param geom Geometry object.
 */
      GridElem(const Lucee::GridGeometry<3, REAL>& geom)
        : FaceElem<REAL>(geom.faceCenter, geom.faceArea)
      {
      }

/**
 * Increment vertex pointer by one element.
 */
      void incr() const { FaceElem<REAL>::incr(); }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
      bool atEnd() const { return FaceElem<REAL>::atEnd(); }
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
      GridElem(const Lucee::GridGeometry<3, REAL>& geom)
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
