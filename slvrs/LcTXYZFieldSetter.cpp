/**
 * @file	LcTXYZFieldSetter.cpp
 *
 * @brief	Linear combiner updater.
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
#include <LcExcept.h>
#include <LcField.h>
#include <LcObjCreator.h>
#include <LcStructuredGridBase.h>
#include <LcTXYZFieldSetter.h>

namespace Lucee
{
// set ids for creators
  template <> const char *TXYZFieldSetter<1>::id = "TXYZFieldSetter1D";
  template <> const char *TXYZFieldSetter<2>::id = "TXYZFieldSetter2D";
  template <> const char *TXYZFieldSetter<3>::id = "TXYZFieldSetter3D";

  template <unsigned NDIM>
  TXYZFieldSetter<NDIM>::TXYZFieldSetter()
    : Lucee::UpdaterIfc(), func(0), isOwner(false)
  {
  }

  template <unsigned NDIM>
  TXYZFieldSetter<NDIM>::~TXYZFieldSetter()
  {
    if (isOwner)
      delete func;
  }

  template <unsigned NDIM>
  void
  TXYZFieldSetter<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    if (tbl.hasTable("func") == false)
      throw Lucee::Except("TXYZFieldSetter<NDIM>::readInput: Must provide a 'func' block");

// create function
    Lucee::LuaTable funcTbl = tbl.getTable("func");
    std::string kind = funcTbl.getKind();
    func = Lucee::ObjCreator<Lucee::FunctionIfc>::getNew(kind);
    func->readInput(funcTbl);

    isOwner = true; // we own it as we made it
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  TXYZFieldSetter<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid = this->template
      getGrid<Lucee::StructuredGridBase<NDIM> >();
// get hold of output array
    Lucee::Field<NDIM, double>& outFld = this->template
      getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> outPtr = outFld.createPtr();

    std::vector<double> txyz(4);
    txyz[0] = t;

    Lucee::RowMajorSequencer<NDIM> seq(outFld.getExtRegion());
    while (seq.step())
    {
// get centroid coordinates
      grid.setIndex(seq.getIndex());
      grid.getCentriod(&txyz[1]); // first txyz element is time

// evaluate function
      std::vector<double> res = func->eval(txyz);

// set output array
      outFld.setPtr(outPtr, seq.getIndex());
      for (unsigned n=0; n<res.size(); ++n)
        outPtr[n] = res[n];
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  TXYZFieldSetter<NDIM>::declareTypes()
  {
// only one output field is expected
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  TXYZFieldSetter<NDIM>::setFunObj(Lucee::FunctionIfc& f)
  {
    if (isOwner)
      delete func;
    func = &f;
    isOwner = false;
  }

// instantiations
  template class TXYZFieldSetter<1>;
  template class TXYZFieldSetter<2>;
  template class TXYZFieldSetter<3>;
}
