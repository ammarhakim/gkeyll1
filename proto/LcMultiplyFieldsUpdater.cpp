/**
 * @file	LcMultiplyFieldsUpdater.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMultiplyFieldsUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *MultiplyFieldsUpdater<1>::id = "MultiplyFields1D";
  template <> const char *MultiplyFieldsUpdater<2>::id = "MultiplyFields2D";
  template <> const char *MultiplyFieldsUpdater<3>::id = "MultiplyFields3D";

  template <unsigned NDIM>
  MultiplyFieldsUpdater<NDIM>::MultiplyFieldsUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  MultiplyFieldsUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
  }

  template <unsigned NDIM>
  void
  MultiplyFieldsUpdater<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  void
  MultiplyFieldsUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  MultiplyFieldsUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    Lucee::Field<NDIM, double>& outFld = this->getOut<Lucee::Field<NDIM, double> >(0);
    Lucee::FieldPtr<double> outFldPtr = outFld.createPtr();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);

    int numFields = this->getNumInpVars();
    int idx[NDIM];

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      outFld.setPtr(outFldPtr, idx);
      outFldPtr = 1.0;
    }

    for (unsigned nf=0; nf<numFields; nf++)
    {
      const Lucee::Field<NDIM, double>& inFld = this->getInp<Lucee::Field<NDIM, double> >(nf);
      Lucee::ConstFieldPtr<double> inFldPtr = inFld.createConstPtr();

// loop, performing integration
      while (seq.step())
      {
	seq.fillWithIndex(idx);
	inFld.setPtr(inFldPtr, idx);
	outFld.setPtr(outFldPtr, idx);
	for (int k=0; k<inFld.getNumComponents(); k++)
	  outFldPtr[k] = outFldPtr[k] * inFldPtr[k];
      }
    }
    return Lucee::UpdaterStatus();
  }

// instantiations
  template class Lucee::MultiplyFieldsUpdater<1>;
  template class Lucee::MultiplyFieldsUpdater<2>;
  template class Lucee::MultiplyFieldsUpdater<3>;
}
