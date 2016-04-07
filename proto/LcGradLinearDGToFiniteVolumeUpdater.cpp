/**
 * @file	LcGradLinearDGToFiniteVolumeUpdater.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcGradLinearDGToFiniteVolumeUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *GradLinearDGToFiniteVolumeUpdater<1>::id = "GradLinearDGToFiniteVolume1D";
  template <> const char *GradLinearDGToFiniteVolumeUpdater<2>::id = "GradLinearDGToFiniteVolume2D";
  template <> const char *GradLinearDGToFiniteVolumeUpdater<3>::id = "GradLinearDGToFiniteVolume3D";

  template <unsigned NDIM>
  GradLinearDGToFiniteVolumeUpdater<NDIM>::GradLinearDGToFiniteVolumeUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  GradLinearDGToFiniteVolumeUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    extrapolateNodes = false;
    if (tbl.hasBool("extrapolateDomainBoundaryNodes"))
      extrapolateNodes = tbl.getBool("extrapolateDomainBoundaryNodes");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  GradLinearDGToFiniteVolumeUpdater<NDIM>::update(double t)
  {
    if (NDIM>1)
      throw Lucee::Except("GradLinearDGToFiniteVolumeUpdater::update: Only works in 1D at present!");
    
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& dgFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& fvFld = this->getOut<Lucee::Field<NDIM, double> >(0);

    double dx1 = 1/grid.getDx(0);
    
    Lucee::ConstFieldPtr<double> dgPtr = dgFld.createConstPtr();    
    Lucee::FieldPtr<double> fvPtr = fvFld.createPtr();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
    {
      fvFld.setPtr(fvPtr, i);
      dgFld.setPtr(dgPtr, i);
      fvPtr[0] = dx1*(dgPtr[1]-dgPtr[0]); // simple central differences
    }
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  GradLinearDGToFiniteVolumeUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class Lucee::GradLinearDGToFiniteVolumeUpdater<1>;
  template class Lucee::GradLinearDGToFiniteVolumeUpdater<2>;
  template class Lucee::GradLinearDGToFiniteVolumeUpdater<3>;
}
