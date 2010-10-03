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
  template <unsigned NDIM, typename T>
  void
  FieldFactory<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
// get name of grid on which field lives
    if (tbl.template hasObject<Lucee::StructuredGridBase<NDIM> >("onGrid"))
      onGrid = &tbl.template
        getObject<Lucee::StructuredGridBase<NDIM> > ("onGrid");
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

  template <unsigned NDIM, typename T>
  Lucee::Field<NDIM, T>*
  FieldFactory<NDIM, T>::create()
  {
// local region for field
    typename Lucee::Region<NDIM, int> localRgn = onGrid->getLocalBox();
// create new field and return pointer to it
    return new Field<NDIM, T>(localRgn, numComponents, lowerGhost, upperGhost);
  }

// instantiations
  template class FieldFactory<1, int>;
  template class FieldFactory<2, int>;
  template class FieldFactory<3, int>;
  template class FieldFactory<4, int>;
  template class FieldFactory<5, int>;
  template class FieldFactory<6, int>;
  template class FieldFactory<7, int>;

  template class FieldFactory<1, float>;
  template class FieldFactory<2, float>;
  template class FieldFactory<3, float>;
  template class FieldFactory<4, float>;
  template class FieldFactory<5, float>;
  template class FieldFactory<6, float>;
  template class FieldFactory<7, float>;

  template class FieldFactory<1, double>;
  template class FieldFactory<2, double>;
  template class FieldFactory<3, double>;
  template class FieldFactory<4, double>;
  template class FieldFactory<5, double>;
  template class FieldFactory<6, double>;
  template class FieldFactory<7, double>;
}
