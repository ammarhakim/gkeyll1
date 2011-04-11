/**
 * @file	LcUnstructGrid.h
 *
 * @brief	Unstructured grid class.
 *
 * @version	$Id$
 */

#ifndef LC_UNSTRUCT_GRID_H
#define LC_UNSTRUCT_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridIfc.h>
#include <LcUnstructConnectivity.h>
#include <LcUnstructGeometry.h>
#include <LcUnstructGridCreator.h>
#include <LcUnstructGridElems.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to represent an unstructured grid.
 */
  template <typename REAL>
  class UnstructGrid : public Lucee::GridIfc
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create an unstructured grid. This grid can not be used unless
 * constructFromCreator() method is called.
 */
      UnstructGrid();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Construct grid from supplied creator.
 *
 * @param ctor Creator to construct grid from.
 */
      void constructFromCreator(const Lucee::UnstructGridCreator<REAL>& ctor);

/**
 * Get number of vertices in grid.
 *
 * @return Number of vertices in grid.
 */
      unsigned getNumVertices() const;

/**
 * Get number of cells in grid. A cell is defined as an element with
 * 'ndim' dimension. I.e. in 3D a cell could be a tetrahedron while in
 * 2D it could be a triangle.
 *
 * @return Number of cells in grid.
 */
      unsigned getNumCells() const;

/**
 * Get number of triangles in grid. This will return 0 for 3D grids,
 * i.e. triangular faces in 3D are not counted as cells.
 *
 * @return number of triangles.
 */
      unsigned getNumTriCells() const;

/**
 * Get number of quadrilaterals in grid.  This will return 0 for 3D
 * grids, i.e. quad faces in 3D are not counted as cells.
 *
 * @return number of quadrilaterals.
 */
      unsigned getNumQuadCells() const;

/**
 * Get number of tetrahedra in grid.
 *
 * @return number of tetrahedra.
 */
      unsigned getNumTetCells() const;

/**
 * Get number of hexahedra in grid.
 *
 * @return number of hexahedra.
 */
      unsigned getNumHexCells() const;

/**
 * Write grid to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the grid as it should appear in output.
 * @return node to which data was written.
 */
      virtual Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
        const std::string& nm);

/**
 * Iterator class to allow iteration over various grid elements. The
 * element dimension is specified as a template parameter.
 */
      template <unsigned NDIM>
      class ElemIterator
      {
// declare unstruct grid as friend so it fiddle with privates
          template <typename R> friend class UnstructGrid;
        public:
/**
 * Create an iterator given grid.
 *
 * @param grid Grid to create iterator from.
 */
          ElemIterator(const UnstructGrid<REAL>& grid)
            : currElem(grid.geometry.vcoords)
          {
          }

/**
 * Return pointer to current element.
 *
 * @return pointer to current element.
 */
          const Lucee::GridElem<REAL, NDIM>* operator->() const
          {
            return &currElem;
          }

/**
 * Increment iterator by single location. (Prefix operator)
 *
 * @return reference to iterator.
 */
          ElemIterator<NDIM>& operator++()
          {
            currElem.incr(); // bump internal location of element
            return *this;
          }

/**
 * Increment iterator by single location. (Postfix operator)
 *
 * @return reference to iterator.
 */
          ElemIterator<NDIM>& operator++(int)
          {
            currElem.incr(); // bump internal location of element
            return *this;
          }

/**
 * Are we at end of iteration?
 *
 * @return true
 */
          bool atEnd() const { return currElem.atEnd(); }
          
        private:
/** Current grid element pointed to by iterator */
          Lucee::GridElem<REAL, NDIM> currElem;
      };

/**
 * Iterator class to allow iteration over various incidence
 * relations. The particular incidence relation is specified as
 * template parameters: Give MDIM and IDIM, this class allows stepping
 * over MDIM->IDIM incidence.
 */
      template <unsigned MDIM, unsigned IDIM>
      class IncidenceIterator : public ElemIterator<MDIM>
      {
// declare unstruct grid as friend so it fiddle with privates
          template <typename R> friend class UnstructGrid;

        public:
/**
 * Create an iterator given grid.
 *
 * @param grid Grid to create iterator from.
 */
          IncidenceIterator(const UnstructGrid<REAL>& grid)
            : ElemIterator<MDIM>(grid)
          {
          }
      };

    private:
/** Dimension of grid */
      unsigned ndim;
/** Geometry information (assume points are in 3d even for 2d mesh) */
      Lucee::UnstructGeometry<3, REAL> geometry;
/** Location 4*d+dprime indicates connectivity d->dprime is stored */
      std::vector<bool> ddprime;
/** Location 4*d+dprime stores d->dprime connectivity information */
      std::vector<Lucee::UnstructConnectivity> connectivity;
/** number of cells of each type */
      std::map<short, unsigned> cellCount;
  };
}

#endif // LC_UNSTRUCT_GRID_H
