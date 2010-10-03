/**
 * @file	LcFieldFactory.cpp
 *
 * @brief	A factory to make fields that live on grids.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcFieldFactory.h>

namespace Lucee
{
  template <unsigned NDIM>
  void
  FieldFactory<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// get name of grid on which field lives
    if (tbl.hasString("onGrid"))
      onGrid = tbl.getString("onGrid");
    else
      throw Lucee::Except("FieldFactory::readInput: must specify 'onGrid', the grid on which field lives");

    numComponents = 1;
// get number of components in field
    if (tbl.hasNumber("numComponents"))
      numComponents = tbl.getNumber("numComponents");

    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
// read number of ghost cells
    if (tbl.hasNumVec("ghost"))
    {
      std::vector<double> gstDbl = tbl.getNumVec("ghost");
      if (gstDbl.size() != 2)
        throw Lucee::Except("FieldFactory::readInput: must specify exactly two entries in 'ghost' table");
      for (unsigned i=0; i<NDIM; ++i)
      {
        lowerGhost[i] = gstDbl[0];
        upperGhost[i] = gstDbl[1];
      }
    }
  }

  template <unsigned NDIM>
  Lucee::Field<NDIM, double>*
  FieldFactory<NDIM>::create()
  {
    throw Lucee::Except("FieldFactory<NDIM>::create: NOT IMPLEMENTED!");
    return 0;
// // cast to solver assembly
//     const Lucee::SolverAssembly& slvrAssembly
//       = dynamic_cast<const Lucee::SolverAssembly&>(solver);
// // get grid
//     const Lucee::StructuredGridBase<NDIM>& grid
//       = slvrAssembly.template getConstGrid<Lucee::StructuredGridBase<NDIM> >(onGrid);
// // local region for field 
//     typename Lucee::Region<NDIM, int> localRgn = grid.getLocalBox();
// // create new field and return pointer to it
//     return new Field<NDIM, double>(localRgn,
//       numComponents, lowerGhost, upperGhost);
  }

// instantiations
  template class FieldFactory<1>;
  template class FieldFactory<2>;
  template class FieldFactory<3>;
}
