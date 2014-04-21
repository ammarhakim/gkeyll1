/**
 * @file LcNonUniRectCartGrid.h
 *
 * @brief A grid allowing for non-uniform spacing, but otherwise rectangular
 */

#ifndef LC_NON_UNI_RECT_CART_GRID_H
#define LC_NON_UNI_RECT_CART_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridGeometry.h>
#include <LcRowMajorIndexer.h>
#include <LcStructuredGridBase.h>
#include <LcVec3.h>

namespace Lucee
{
/**
 * Represents a non-uniform but rectangular Cartesian grid of
 * arbitrary dimensions.
 */
  template <unsigned NDIM>
  class NonUniRectCartGrid : public Lucee::StructuredGridBase<NDIM>
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Default ctor
 */
      NonUniRectCartGrid();

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
      virtual void getCentroid(double xc[3]) const;

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
      virtual TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

    private:
  };
}

#endif // LC_NON_UNI_RECT_CART_GRID_H
