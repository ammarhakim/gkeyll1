/**
 * @file	LcRecordFieldDerivInCell.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcRecordFieldDerivInCell.h>
#include <LcStructuredGridBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *RecordFieldDerivInCell<1>::id = "RecordFieldDerivInCell1D";
  template <> const char *RecordFieldDerivInCell<2>::id = "RecordFieldDerivInCell2D";
  template <> const char *RecordFieldDerivInCell<3>::id = "RecordFieldDerivInCell3D";

  template <unsigned NDIM>
  RecordFieldDerivInCell<NDIM>::RecordFieldDerivInCell()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  RecordFieldDerivInCell<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of cell to record data
    if (tbl.hasNumVec("cellIndex"))
    {
      std::vector<double> myIdx = tbl.getNumVec("cellIndex");
      if (myIdx.size() != NDIM)
      {
        Lucee::Except lce("Lucee::RecordFieldDerivInCell: Must specify 'cellIndex' with ");
        lce << NDIM << " entries. Provided " << myIdx.size() <<  " entries instead";
        throw lce;
      }
      for (unsigned i=0; i<NDIM; ++i)
        cellIdx[i] = (int) myIdx[i];
    }
    else
    {
      Lucee::Except lce("Lucee::RecordFieldDerivInCell: Must specify 'cellIndex' ");
      throw lce;
    }

// get hold of derivative to compute
    if (tbl.hasNumber("order"))
      order = tbl.hasNumber("order");
    else
    {
      Lucee::Except lce("Lucee::RecordFieldDerivInCell: Must specify 'order' ");
      throw lce;
    }
    if (order > 2)
    {
      Lucee::Except lce("Lucee::RecordFieldDerivInCell: Only 1st od 2nd order derivatives are supported.");
      throw lce;
    }
  }

  template <unsigned NDIM>
  void
  RecordFieldDerivInCell<NDIM>::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  RecordFieldDerivInCell<NDIM>::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
// get input array
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
// get output dynVector
    Lucee::DynVector<double>& fldVal = this->getOut<Lucee::DynVector<double> >(0);

    unsigned nc = fld.getNumComponents();
    int cellIdx_p[NDIM], cellIdx_m[NDIM];

    std::vector<double> localVals(NDIM*nc);
    for (unsigned i=0; i<localVals.size(); ++i)
      localVals[i] = 0.0;

// local region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
// record field value if inside domain
    if (localRgn.isInside(cellIdx))
    {
      Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_p = fld.createConstPtr();
      Lucee::ConstFieldPtr<double> fldPtr_m = fld.createConstPtr();

      fld.setPtr(fldPtr, cellIdx);
// compute derivatives in each direction
      for (unsigned d=0; d<NDIM; ++d)
      {
        for (unsigned dd=0; dd<NDIM; ++dd)
        {
          cellIdx_p[dd] = cellIdx_m[dd] = cellIdx[dd];
        }
        cellIdx_p[d] += 1;
        cellIdx_m[d] -= 1;

        fld.setPtr(fldPtr_p, cellIdx_p);
        fld.setPtr(fldPtr_m, cellIdx_m);

        double dx = grid.getDx(d);
        if (order == 1)
          for (unsigned i=0; i<nc; ++i)
            localVals[d*nc+i] = (fldPtr_p[i]-fldPtr_m[i])/(2*dx);
        else if (order == 2)
          for (unsigned i=0; i<nc; ++i)
            localVals[d*nc+i] = (fldPtr_p[i]-2*fldPtr[i]+fldPtr_m[i])/(dx*dx);
      }
    }

// filled up from all reduce across all processors
    std::vector<double> commLocalVals(NDIM*fld.getNumComponents());
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
// sum across all processors (this works because except for one
// processor, all others have zeros)
    comm->allreduce(NDIM*fld.getNumComponents(), localVals, commLocalVals, TX_SUM);

// push value into dynvector
    fldVal.appendData(t, commLocalVals);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  RecordFieldDerivInCell<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

// instantiations
  template class Lucee::RecordFieldDerivInCell<1>;
  template class Lucee::RecordFieldDerivInCell<2>;
  template class Lucee::RecordFieldDerivInCell<3>;
}
