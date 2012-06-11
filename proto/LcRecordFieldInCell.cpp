/**
 * @file	LcRecordFieldInCell.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcRecordFieldInCell.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *RecordFieldInCell<1>::id = "RecordFieldInCell1D";
  template <> const char *RecordFieldInCell<2>::id = "RecordFieldInCell2D";
  template <> const char *RecordFieldInCell<3>::id = "RecordFieldInCell3D";

  template <unsigned NDIM>
  RecordFieldInCell<NDIM>::RecordFieldInCell()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RecordFieldInCell<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of cell to record data
    if (tbl.hasNumVec("cellIndex"))
    {
      std::vector<double> myIdx = tbl.getNumVec("cellIndex");
      if (myIdx.size() != NDIM)
      {
        Lucee::Except lce("Lucee::RecordFieldInCell: Must specify 'cellIndex' with ");
        lce << NDIM << " entries. Provided " << myIdx.size() <<  " entries instead";
        throw lce;
      }
      for (unsigned i=0; i<NDIM; ++i)
        cellIdx[i] = (int) myIdx[i];
    }
    else
    {
      Lucee::Except lce("Lucee::RecordFieldInCell: Must specify 'cellIndex' ");
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  RecordFieldInCell<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RecordFieldInCell<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input array
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
// get output dynVector
    Lucee::DynVector<double>& fldVal = this->getOut<Lucee::DynVector<double> >(0);

    std::vector<double> localVals(fld.getNumComponents());
    for (unsigned i=0; i<localVals.size(); ++i)
      localVals[i] = 0.0;

// local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
// record field value if inside domain
    if (localRgn.isInside(cellIdx))
    {
      Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
      fld.setPtr(fldPtr, cellIdx);
      for (unsigned i=0; i<localVals.size(); ++i)
        localVals[i] = fldPtr[i];
    }

// TODO TODO all reduce so everyone has data: should an all-reduce be
// done before dump or here?

// push value into dynvector
    fldVal.appendData(t, localVals);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RecordFieldInCell<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

// instantiations
  template class Lucee::RecordFieldInCell<1>;
  template class Lucee::RecordFieldInCell<2>;
  template class Lucee::RecordFieldInCell<3>;
}
