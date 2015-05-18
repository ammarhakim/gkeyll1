/**
 * @file	LcUnstructuredGrid.h
 *
 * @brief	Base class for moab unstructred grids
 */

#ifndef LC_UNSTRUCTURED_GRID_H
#define LC_UNSTRUCTURED_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDecompRegionCalcIfc.h>
#include <LcGridIfc.h>
#include <LcRegion.h>

// boost includes
#include <boost/shared_ptr.hpp>

namespace Lucee
{
/**
 * A base class to represent a single-block rectangular body-fitted
 * grid. The indexing scheme is designed such that the cell owns the
 * faces on the "left", "bottom" and "back" (0, 1, 2
 * directions). I.e. the cell index also can be used to get surface
 * variables on these faces using the proper direction index.
 */
  template <unsigned NDIM>
  class UnstructuredGrid : public Lucee::GridIfc
  {
    public:
/**
 * Destroy object.
 */
      virtual ~UnstructuredGrid();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Get number of cells in given direction.
 *
 * @param dir Direction in which number of cells is needed.
 * @return Number of cells.
 */
      unsigned getNumCells(unsigned dir) const;


/**
 * Set the current cell location in grid to (i).
 *
 * @param i Index location into grid.
 */
      void setIndex(int i) const;

/**
 * Return coordinates in physical space of cell centroid.
 *
 * @param xc On output, centroid of cell.
 */
      virtual void getCentroid(double xc[]) const = 0;

/**
 * Return coordinates in physical space of bottom left node.
 *
 * @param xc On output, vertex coordinate of cell.
 */
      virtual void getVertex(double xc[]) const = 0;

/**
 * Return volume of cell.
 *
 * @return Cell volume.
 */
      virtual double getVolume() const = 0;

/**
 * Return physical surface area of face perpendicular (in
 * computational space) to specified direction.
 *
 * @param dir Direction perpendicular to face.
 * @return surface area of face.
 */
      virtual double getSurfArea(unsigned dir) const = 0;

/**
 * Get the coordinate system attached to the surface perpendicular (in
 * computational space) to specified direction. The coordinate system
 * is defined by three unit vectors, norm, tan1 and tan2. The norm
 * vector is normal to the surface and points into the cell. The
 * vectors tan1 and tan2 lie in the plane of the surface. The system
 * is such that tan1 x tan2 = norm.
 *
 * @param dir Direction perpendicular to face.
 * @param norm On output, normal to face.
 * @param tan1 On output, first tangent to face.
 * @param tan2 On output, second tangent to face.
 * 
 */
      virtual void getSurfCoordSys(unsigned dir, double norm[],
        double tan1[], double tan2[]) const = 0;

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);


/**
 * Lua callable method to get global shape
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetShape(lua_State *L);


    protected:
/**
 * Default ctor: only derived classes can make default objects.
 */
      UnstructuredGrid();

/**
 * Create a structured grid on specified region.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param globalRgn Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      UnstructuredGrid(const Lucee::Region<NDIM, int>& globalRgn,
        const Lucee::Region<NDIM, double>& compSpace);

/**
 * Create a structured grid using specified decomposed region.
 *
 * @param dcmpRgn Decomposed region to use.
 * @param compSpace Region in computation space.
 */
      UnstructuredGrid(const Lucee::DecompRegion<NDIM>& dcmpRgn,
        const Lucee::Region<NDIM, double>& compSpace);

/**
 * Set structured grid from supplied one.
 *
 * @param sg Structured grid to assign from.
 * @return reference to this object.
 */
     UnstructuredGrid<NDIM>& operator=(const UnstructuredGrid<NDIM>& sg);


/** Index into current cell */
      mutable int currIdx[NDIM];

    private:

      //Not sure if these belong here or in geometry.
      moab::Interface* mb;

      moab::Range faces;
      moab::Range cells;
      moab::Range vertices;
      moab::Range edges;

  };
}

#endif //  LC_STRUCTURED_GRID_BASE_H
