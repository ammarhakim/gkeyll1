/**
 * @file	LcOverlappingFieldSplit.cpp
 *
 * @brief	Updater to copy data from global domain to two split domains
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcOverlappingFieldSplit.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <> const char *OverlappingFieldSplit<1>::id = "OverlappingFieldSplit1D";
  template <> const char *OverlappingFieldSplit<2>::id = "OverlappingFieldSplit2D";
  template <> const char *OverlappingFieldSplit<3>::id = "OverlappingFieldSplit3D";

  template <unsigned NDIM>
  OverlappingFieldSplit<NDIM>::OverlappingFieldSplit()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  OverlappingFieldSplit<NDIM>::~OverlappingFieldSplit()
  {
  }

  template <unsigned NDIM>
  void
  OverlappingFieldSplit<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method first
    Lucee::UpdaterIfc::readInput(tbl);

    numCells = tbl.getNumber("numOverlappingCells");
    dir = tbl.getNumber("dir");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  OverlappingFieldSplit<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    const Lucee::Field<NDIM, double>& qFull = this->getInp<Lucee::Field<NDIM, double> >(0);
    
    Lucee::Field<NDIM, double>& qLeft = this->getOut<Lucee::Field<NDIM, double> >(0);;
    Lucee::Field<NDIM, double>& qRight = this->getOut<Lucee::Field<NDIM, double> >(1);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  OverlappingFieldSplit<NDIM>::declareTypes()
  {
// expect one input field
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
// expect two output field
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));    
  }

// instantiations
  template class OverlappingFieldSplit<1>;
  template class OverlappingFieldSplit<2>;
  template class OverlappingFieldSplit<3>;
}
