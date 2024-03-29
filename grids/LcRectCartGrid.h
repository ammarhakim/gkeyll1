/**
 * @file	LcRectCartGrid.h
 *
 * @brief	A rectangular cartesian grid.
 */

#ifndef LC_RECT_CART_GRID_H
#define LC_RECT_CART_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/**
 * Represents a rectangular cartesian grid of arbitrary dimensions.
 */
  template <unsigned NDIM>
  class RectCartGrid : public Lucee::StructuredGridBase<NDIM>
  {
// Number of components for coordinate arrays etc.
      static const unsigned NC = Lucee::CDIM<NDIM>::N;
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Default ctor
 */
      RectCartGrid();

/**
 * Create a new rectangular Cartesian grid on specified region.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param globalBox Global index region for this grid.
 * @param physBox Rectangular region in space.
 */
      RectCartGrid(const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& physBox);

/**
 * Create a new rectangular grid using specified decomposed region.
 *
 * @param dcmpRgn Decomposed region to use.
 * @param physBox Rectangular region in space.
 */
      RectCartGrid(const Lucee::DecompRegion<NDIM>& dcmpRgn,
        const Lucee::Region<NDIM, double>& physSpace);

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Return coordinates in physical space of cell centroid. The
 * setIndex() method must be called before this to set the current
 * cell index.
 *
 * @param xc On output, centroid of cell.
 */
      virtual void getCentroid(double xc[]) const;

/**
 * Return coordinates in physical space of bottom left node. The
 * setIndex() method must be called before this to set the current
 * cell index.
 *
 * @param xc On output, vertex coordinate of cell.
 */
      virtual void getVertex(double xc[]) const;

/**
 * Return volume of cell. The setIndex() method must be called before
 * this to set the current cell index.
 *
 * @return Cell volume.
 */
      virtual double getVolume() const;

/**
 * Return physical surface area of face perpendicular (in
 * computational space) to specified direction. The setIndex() method
 * must be called before this to set the current cell index.
 *
 * @param dir Direction perpendicular to face.
 * @return surface area of face.
 */
      virtual double getSurfArea(unsigned dir) const;

/**
 * Get the coordinate system attached to the surface perpendicular (in
 * computational space) to specified direction. The coordinate system
 * is defined by three unit vectors, norm, tan1 and tan2. The norm
 * vector is normal to the surface and points into the cell. The
 * vectors tan1 and tan2 lie in the plane of the surface. The system
 * is such that tan1 x tan2 = norm. The setIndex() method must be
 * called before this to set the current cell index.
 *
 * @param dir Direction perpendicular to face.
 * @param norm On output, normal to face.
 * @param tan1 On output, first tangent to face.
 * @param tan2 On output, second tangent to face.
 * 
 */
      virtual void getSurfCoordSys(unsigned dir, double norm[],
        double tan1[], double tan2[]) const;

/**
 * Write grid to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the grid as it should appear in output.
 * @return node to which data was written.
 */
      virtual TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
/** Grid spacing in each direction */
      double dx[Lucee::CDIM<NDIM>::N];
/** Volume of each cell */
      double cellVolume;
/** Number of cells in domain */
      std::vector<double> cells;
/** Lower coordinates of space */
      std::vector<double> lower;
/** Upper coordinates of space */
      std::vector<double> upper;

/**
 * Copy from supplied rectangular grid.
 *
 * @param rg Rectangular grid to copy from.
 * @return reference to this object.
 */
      RectCartGrid<NDIM>& operator=(const RectCartGrid<NDIM>& rg);

/**
 * Compute basic geometric information needed by grid.
 */
      void calcGeometry();
  };
}

#endif // LC_RECT_CART_GRID_H
