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
 * Write grid to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the grid as it should appear in output.
 * @return node to which data was written.
 */
      virtual Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
        const std::string& nm);

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
