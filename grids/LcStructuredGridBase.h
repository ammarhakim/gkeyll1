/**
 * @file	LcStructuredGridBase.h
 *
 * @brief	Base class for body fitted grid in arbitrary dimensions.
 */

#ifndef LC_STRUCTURED_GRID_BASE_H
#define LC_STRUCTURED_GRID_BASE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridIfc.h>
#include <LcRegion.h>

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
  class StructuredGridBase : public Lucee::GridIfc
  {
    public:
/**
 * Destroy object.
 */
      virtual ~StructuredGridBase();

/**
 * Get number of cells in given direction.
 *
 * @param dir Direction in which number of cells is needed.
 * @return Number of cells.
 */
      unsigned getNumCells(unsigned dir) const;

/**
 * Return global region indexed by grid.
 *
 * @return global region indexed by grid.
 */
      Lucee::Region<NDIM, int> getGlobalBox() const;

/**
 * Return local region indexed by grid.
 *
 * @return  local region indexed by grid.
 */
      Lucee::Region<NDIM, int> getLocalBox() const;

/**
 * Return region in computational space
 *
 * @return region in computational space
 */
      Lucee::Region<NDIM, double> getComputationalSpace() const;

/**
 * Return cell spacing in computational space in specifed direction.
 *
 * @param dir Direction in which spacing is needed.
 * @return spacing in computational space in direction 'dir'
 */
      double getDx(unsigned dir) const;

/**
 * Set the current cell location in grid to (i).
 *
 * @param i Index location into grid.
 */
      void setIndex(int i) const;

/**
 * Set the current cell location in grid to (i, j).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 */
      void setIndex(int i, int j) const;

/**
 * Set the current cell location in grid to (i, j, k).
 *
 * @param i Index location into grid.
 * @param j Index location into grid.
 * @param k Index location into grid.
 */
      void setIndex(int i, int j, int k) const;

/**
 * Set the current cell location in grid to specified index.
 *
 * @param idx Index location into grid.
 */
      void setIndex(const int idx[NDIM]) const;

/**
 * Return coordinates in physical space of cell centroid.
 *
 * @param xc On output, centroid of cell.
 */
      virtual void getCentroid(double xc[3]) const = 0;

/**
 * Return coordinates in physical space of bottom left node.
 *
 * @param xc On output, vertex coordinate of cell.
 */
      virtual void getVertex(double xc[3]) const = 0;

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
 * computational space) to specified direction. The corrdinate system
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
      virtual void getSurfCoordSys(unsigned dir, double norm[3],
        double tan1[3], double tan2[3]) const = 0;

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method to get lower index of local box.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetLocalLower(lua_State *L);

/**
 * Lua callable method to get upper index of local box.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetLocalUpper(lua_State *L);

/**
 * Lua callable method to get lower index of global box.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetGlobalLower(lua_State *L);

/**
 * Lua callable method to get upper index of global box.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetGlobalUpper(lua_State *L);

/**
 * Lua callable method to get coordinate of lower corner of domain in
 * computational space.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetLowerCoord(lua_State *L);

/**
 * Lua callable method to get coordinate of upper corner of domain in
 * computational space.
 *
 * @param L Lua state to work with.
 * @return number of return values.
 */
      static int luaGetUpperCoord(lua_State *L);

    protected:
/**
 * Default ctor: only derived classes can make default objects.
 */
      StructuredGridBase();

/**
 * Create a new body-fitted grid on specified region. In serial the
 * local and global boxes coincide. In parallel, the globalBox
 * represents the full grid, while the localBox represents the portion
 * of the grid handled by the rank the grid lives on. The full
 * computational domain (i.e. indexed by the global region) is
 * specified using the Region object compSpace.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param localBox Local index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      StructuredGridBase(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& compSpace);

/**
 * Set structured grid from supplied one.
 *
 * @param sg Structured grid to assign from.
 * @return reference to this object.
 */
      StructuredGridBase<NDIM>& operator=(const StructuredGridBase<NDIM>& sg);

/**
 * Set data required to build grid.
 *
 * @param localBox Local index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param compSpace Region in computation space.
 */
      void setGridData(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& globalBox, const Lucee::Region<NDIM, double>& compSpace);

/** Index into current cell */
      mutable int currIdx[NDIM];
/** Local region indexed by grid */
      Lucee::Region<NDIM, int> localBox;
/** Global region indexed by grid */
      Lucee::Region<NDIM, int> globalBox;
/** Global region spanned by grid in computational space */
      Lucee::Region<NDIM, double> compSpace;
  };
}

#endif //  LC_STRUCTURED_GRID_BASE_H
