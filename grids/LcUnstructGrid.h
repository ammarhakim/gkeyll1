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
        public:
/**
 * Return pointer to current element.
 *
 * @return pointer to current element.
 */
          const Lucee::GridElem<REAL, NDIM>* operator->() const
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
            currElem->incr(); // bump internal location of element
            return *this;
          }

/**
 * Increment iterator by single location. (Postfix operator)
 *
 * @return reference to iterator.
 */
          ElemIterator<NDIM>& operator++(int)
          {
            currElem->incr(); // bump internal location of element
            return *this;
          }
          
        private:
/** Current grid element pointed to by iterator */
          Lucee::GridElem<REAL, NDIM> *currElem;
      };

    private:
/** Dimension of grid */
      unsigned ndim;
/** Geometry information (assume 3 nodal coordinates) */
      Lucee::UnstructGeometry<3, REAL> geometry;
/** Location 4*d+dprime indicates connectivity d->dprime is stored */
      std::vector<bool> ddprime;
/** Location 4*d+dprime stores d->dprime connectivity information */
      std::vector<Lucee::UnstructConnectivity> connectivity;
  };
}

#endif // LC_UNSTRUCT_GRID_H
