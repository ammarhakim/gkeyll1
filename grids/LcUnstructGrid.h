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
 * Get specified grid element.
 *
 * @param idx Index of requested grid element.
 * @return grid element.
 */
      template <unsigned NDIM>
      Lucee::GridElem<REAL, NDIM> getElement(unsigned idx) const
      {
        return Lucee::GridElem<REAL, NDIM>(this->geometry);
      }

/**
 * Iterator class to allow iteration over various grid elements. The
 * element dimension is specified as a template parameter.
 */
      template <unsigned NDIM>
      class ElemIterator
      {
        public:
/**
 * Create an iterator given grid.
 *
 * @param grid Grid to create iterator from.
 */
          ElemIterator(const UnstructGrid<REAL>& grid)
            : currElem(grid.geometry)
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
 * Return current grid element.
 *
 * @return Current element.
 */
          Lucee::GridElem<REAL, NDIM> operator*() const
          {
            return currElem;
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
 * template parameters: Given MDIM and IDIM, this class allows
 * stepping over MDIM->IDIM incidence.
 */
      template <unsigned MDIM, unsigned IDIM>
      class IncidenceIterator : public ElemIterator<MDIM>
      {
        public:
/**
 * Create an iterator given grid.
 *
 * @param grid Grid to create iterator from.
 */
          IncidenceIterator(const UnstructGrid<REAL>& grid)
            : ElemIterator<MDIM>(grid), conn(grid.getConnectivity(MDIM,IDIM)), curr(0)
          {
            currIndices = &conn.indices[conn.offsets[curr]];
          }

/**
 * Increment iterator by single location. (Prefix operator)
 *
 * @return reference to iterator.
 */
          IncidenceIterator<MDIM, IDIM>& operator++()
          {
            ElemIterator<MDIM>::operator++();
            curr++;
            currIndices = &conn.indices[conn.offsets[curr]];
            return *this;
          }

/**
 * Increment iterator by single location. (Postfix operator)
 *
 * @return reference to iterator.
 */
          IncidenceIterator<MDIM, IDIM>& operator++(int)
          {
            ElemIterator<MDIM>::operator++(0);
            curr++;
            currIndices = &conn.indices[conn.offsets[curr]];
            return *this;
          }

/**
 * Get current element number.
 *
 * @return current element number.
 */
          unsigned getCurrIndex() const { return curr; }

/**
 * Get number of connections for current element.
 *
 * @return number of connections for current element.
 */
          unsigned getNumConnections() const
          {
            return conn.offsets[curr+1]-conn.offsets[curr];
          }

/**
 * Get index for 'idx' connection of current element.
 *
 * @param idx Index of connection number.
 * @return Index of connected elemented.
 */
          int getIndex(unsigned idx) const
          {
            return currIndices[idx];
          }

        private:
/** Reference to connectivity information */
          const Lucee::UnstructConnectivity& conn;
/** Current element being indexed */
          unsigned curr;
/** Pointer to connectivities for current element */
          const int *currIndices;
      };

    private:
/** Dimension of grid */
      unsigned ndim;
/** Geometry information (assume points are in 3d even for 2d mesh) */
      Lucee::UnstructGeometry<3, REAL> geometry;
/** Location 4*d+dprime indicates connectivity d->dprime is stored */
      std::vector<bool> ddprime;
/** Location 4*d+dprime stores d->dprime connectivity information */
      mutable std::vector<Lucee::UnstructConnectivity> connectivity;
/** number of cells of each type */
      std::map<short, unsigned> cellCount;
/** Cell type */
      std::vector<short> cellType;

/**
 * Method to return reference to specified connectivity,
 * d->dprime. This method will compute the connectivity if not done
 * already.
 *
 * @param d Dimension to connect (d->dprime).
 * @param dprime Dimension to connect to (d->dprime).
 * @return reference to connectivity.
 */
      const Lucee::UnstructConnectivity& getConnectivity(unsigned d, unsigned dprime) const;

/**
 * Compute geometry of cells. This computes area and centroid of each
 * cell in a 2D grid.
 */
      void calcCellGeometry2d();
  };
}

#endif // LC_UNSTRUCT_GRID_H
