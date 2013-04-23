/**
 * @file	LcIntegrateFieldAlongLine.cpp
 *
 * @brief	Updater to integrate a field along a specified line.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcIntegrateFieldAlongLine.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{
  template <> const char *IntegrateFieldAlongLine<1>::id = "IntegrateFieldAlongLine1D";
  template <> const char *IntegrateFieldAlongLine<2>::id = "IntegrateFieldAlongLine2D";
  template <> const char *IntegrateFieldAlongLine<3>::id = "IntegrateFieldAlongLine3D";

  template <unsigned NDIM>
  IntegrateFieldAlongLine<NDIM>::IntegrateFieldAlongLine()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  IntegrateFieldAlongLine<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

// start cell for integration
    std::vector<double> sc = tbl.getNumVec("startCell");
    if (sc.size() != NDIM)
    {
      Lucee::Except lce("IntegrateFieldAlongLine::readInput: 'startCell' must have '");
      lce << NDIM << " entries. Provided " << sc.size() << " instead";
      throw lce;
    }
    for (unsigned i=0; i<NDIM; ++i)
      startCell[i] = (int) sc[i];

// direction to integrate along
    dir = tbl.getNumber("dir");
    if (dir > NDIM-1)
    {
      Lucee::Except lce("IntegrateFieldAlongLine::readInput: 'dir' must be less than ");
      lce << NDIM << ". Provided  " << dir << " instead";
      throw lce;
    }

// number of cells to integrate along
    numCells = (unsigned) tbl.getNumber("numCells");
  }

  template <unsigned NDIM>
  void
  IntegrateFieldAlongLine<NDIM>::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  IntegrateFieldAlongLine<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

// create box for integration
    int upper[NDIM];
    for (unsigned d=0; d<NDIM; ++d)
      upper[NDIM] = startCell[NDIM]+1;
    upper[dir] = startCell[dir]+numCells;
    Lucee::Region<NDIM, int> intRgn(startCell, upper);
// local grid region
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();    
    Lucee::RowMajorSequencer<NDIM> seq(localRgn.intersect(intRgn));

// grid spacing
    double dx = grid.getDx(dir);

    Lucee::ConstFieldPtr<double> fldPtr = fld.createConstPtr();
    int idx[NDIM];
    double localInt = 0.0;
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fld.setPtr(fldPtr, idx);
      localInt += dx*fldPtr[0];
    }

    double volInt = localInt;
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;
    comm->allreduce(1, &localInt, &volInt, TX_SUM);

    std::vector<double> data(1);
    data[0] = volInt;
    fldInt.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  IntegrateFieldAlongLine<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

// instantiations
  template class Lucee::IntegrateFieldAlongLine<1>;
  template class Lucee::IntegrateFieldAlongLine<2>;
  template class Lucee::IntegrateFieldAlongLine<3>;
}

