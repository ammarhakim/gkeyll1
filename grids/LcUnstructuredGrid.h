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

#include <MBCore.hpp>
#include <MBReaderWriterSet.hpp>
#include <MBReadUtilIface.hpp>

namespace Lucee
{

  template<unsigned NDIM>
  class UnstructuredGrid: public Lucee::GridIfc
  {
    public:

      static const char *id;

      /**
       * Default ctor: only derived classes can make default objects.
       */
      UnstructuredGrid();

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
       * Write grid to specified HDF5 file.
       *
       * @param nm Name of file to write.
       */
            virtual void write(const std::string& nm);
#if 1
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
#endif

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
      virtual void getCentroid(double xc[]) const
      {
      }
      ;

      /**
       * Return coordinates in physical space of bottom left node.
       *
       * @param xc On output, vertex coordinate of cell.
       */
      virtual void getVertex(double xc[]) const
      {
      }
      ;

      /**
       * Return volume of cell.
       *
       * @return Cell volume.
       */
      virtual double getVolume() const
      {
        return 0.0;
      }
      ;

      /**
       * Return physical surface area of face perpendicular (in
       * computational space) to specified direction.
       *
       * @param dir Direction perpendicular to face.
       * @return surface area of face.
       */
      virtual double getSurfArea(unsigned dir) const
      {
        return 0.0;
      }
      ;

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
      virtual void getSurfCoordSys(unsigned dir, double norm[], double tan1[],
          double tan2[]) const
      {
      }
      ;

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
       * Set structured grid from supplied one.
       *
       * @param sg Structured grid to assign from.
       * @return reference to this object.
       */
      UnstructuredGrid<NDIM>& operator=(const UnstructuredGrid<NDIM>& sg);

      /** Index into current cell */
      mutable int currIdx[NDIM];

    private:

      std::string filename;
      std::string outFilename;
      std::string writeExtension;

      moab::EntityHandle set;
      char* readOpts;
      //moab::FileOptions writeOpts;
      moab::Interface* mb;

      moab::Range faces;
      moab::Range cells;
      moab::Range vertices;
      moab::Range edges;

  };
}

#endif //  LC_STRUCTURED_GRID_BASE_H
