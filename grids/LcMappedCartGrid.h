/**
 * @file LcMappedCartGrid.h
 *
 * @brief A mapped rectangular grid
 */

#ifndef LC_MAPPED_CART_GRID_H
#define LC_MAPPED_CART_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridGeometry.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
/**
 * Represents a mapped cartesian grid of arbitrary
 * dimensions. I.e. grid is rectangular in computational space.
 */
  template <unsigned NDIM>
  class MappedCartGrid : public Lucee::StructuredGridBase<NDIM>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Default ctor
 */
      MappedCartGrid();

/**
 * Create a new rectangular Cartesian grid on specified region. In
 * serial the local and global boxes coincide. In parallel, the
 * globalBox represents the full grid, while the localBox represents
 * the portion of the grid handled by the rank the grid lives on.
 *
 * To get data from the grid, first set the index into the grid by
 * using the setIndex() method. Then access the needed data for that
 * cell. For structured grids, the edges on the lower side (left,
 * bottom, back) are labeled by the cell index. Further, the lower
 * left corner of each cell is also labelled by the cell index.
 *
 * @param localBox Local index region for this grid.
 * @param localExtBox Local extended index region for this grid.
 * @param globalBox Global index region for this grid.
 * @param compSpace Rectangular region in computational space.
 */
      MappedCartGrid(const Lucee::Region<NDIM, int>& localBox,
        const Lucee::Region<NDIM, int>& localExtBox,
        const Lucee::Region<NDIM, int>& globalBox,
        const Lucee::Region<NDIM, double>& compSpace);

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
      virtual void getCentriod(double xc[3]) const;

/**
 * Return coordinates in physical space of bottom left node. The
 * setIndex() method must be called before this to set the current
 * cell index.
 *
 * @param xc On output, vertex coordinate of cell.
 */
      virtual void getVertex(double xc[3]) const;

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
      virtual void getSurfCoordSys(unsigned dir, double norm[3],
        double tan1[3], double tan2[3]) const;

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
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
/**
 * Copy from supplied rectangular grid.
 *
 * @param rg Rectangular grid to copy from.
 * @return reference to this object.
 */
      MappedCartGrid<NDIM>& operator=(const MappedCartGrid<NDIM>& rg);

/** Computational grid spacing in each direction */
      double dx[3];
/** Volume of computational cell (constant as cells are
 * rectangular) */
      double cellVolume;
/** Extended local region indexed by grid (needed to index the ghost
 * region geometry) */
      Lucee::Region<NDIM, int> localExtBox;
/** Geometry object to store grid geometry */
      Lucee::GridGeometry<NDIM, double> geometry;
  };
}

#endif // LC_MAPPED_CART_GRID_H
