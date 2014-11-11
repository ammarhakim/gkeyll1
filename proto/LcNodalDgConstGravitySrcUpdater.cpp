/**
 * @file	LcNodalDgConstGravitySrcUpdater.cpp
 *
 * @brief	Updater for gravitational source
 */

// lucee includes
#include <LcNodalDgConstGravitySrcUpdater.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

namespace Lucee
{
// set ids for module system
  template <> const char *NodalDgConstGravitySrcUpdater<1>::id = "NodalDgConstGravitySrc1D";
  template <> const char *NodalDgConstGravitySrcUpdater<2>::id = "NodalDgConstGravitySrc2D";
  template <> const char *NodalDgConstGravitySrcUpdater<3>::id = "NodalDgConstGravitySrc3D";

  template <unsigned NDIM>
  void
  NodalDgConstGravitySrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    dir = (unsigned) tbl.getNumber("dir");
    gravity = tbl.getNumber("gravity");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContSrcIncrUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalDgConstGravitySrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  NodalDgConstGravitySrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    double dt = t-this->getCurrTime();
// get fluid
    Lucee::Field<2, double>& fluid = this->getOut<Lucee::Field<2, double> >(0);
    Lucee::FieldPtr<double> fPtr = fluid.createPtr();

    unsigned nlocal = nodalBasis->getNumNodes();

    int dirIdx = dir+1; // offset as first index is rho

    int idx[2];
    Lucee::Region<2, int> localRgn = fluid.getRegion();
    Lucee::RowMajorSequencer<2> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fluid.setPtr(fPtr, idx);

// updated solution at each node
      for (unsigned n=0; n<nlocal; ++n)
      {
        unsigned sIdx = 5*n; // 5 equations

        double rho = fPtr[sIdx+0];
// old contribution to kinetic energy from momentum in 'dir' direction
        double keold = 0.5*fPtr[sIdx+dirIdx]*fPtr[sIdx+dirIdx]/rho;
// update momentum
        fPtr[sIdx+dirIdx] += gravity*rho*dt;
// now update energy to account for updated momentum
        fPtr[sIdx+4] = fPtr[sIdx+4] - keold + 0.5*fPtr[sIdx+dirIdx]*fPtr[sIdx+dirIdx]/rho;
      }
    }
    
    return Lucee::UpdaterStatus();
  }
  
  template <unsigned NDIM>
  void
  NodalDgConstGravitySrcUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

// instantiations
  template class NodalDgConstGravitySrcUpdater<1>;
  template class NodalDgConstGravitySrcUpdater<2>;
  template class NodalDgConstGravitySrcUpdater<3>;
}

