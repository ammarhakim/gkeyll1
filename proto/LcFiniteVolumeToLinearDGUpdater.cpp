/**
 * @file	LcFiniteVolumeToLinearDGUpdater.cpp
 *
 * @brief	Updater to integrate nodal DG/CG field over complete domain
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcFiniteVolumeToLinearDGUpdater.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <> const char *FiniteVolumeToLinearDGUpdater<1>::id = "FiniteVolumeToLinearDG1D";
  template <> const char *FiniteVolumeToLinearDGUpdater<2>::id = "FiniteVolumeToLinearDG2D";
  template <> const char *FiniteVolumeToLinearDGUpdater<3>::id = "FiniteVolumeToLinearDG3D";

  template <unsigned NDIM>
  FiniteVolumeToLinearDGUpdater<NDIM>::FiniteVolumeToLinearDGUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  FiniteVolumeToLinearDGUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    extrapolateNodes = false;
    if (tbl.hasBool("extrapolateDomainBoundaryNodes"))
      extrapolateNodes = tbl.getBool("extrapolateDomainBoundaryNodes");
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  FiniteVolumeToLinearDGUpdater<NDIM>::update(double t)
  {
    if (NDIM>1)
      throw Lucee::Except("FiniteVolumeToLinearDGUpdater::update: Only works in 1D at present!");
    
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    const Lucee::Field<NDIM, double>& fvFld = this->getInp<Lucee::Field<NDIM, double> >(0);
    Lucee::Field<NDIM, double>& dgFld = this->getOut<Lucee::Field<NDIM, double> >(0);

    Lucee::ConstFieldPtr<double> fvPtr = fvFld.createConstPtr();
    Lucee::ConstFieldPtr<double> fvPtr_m = fvFld.createConstPtr();
    Lucee::ConstFieldPtr<double> fvPtr_p = fvFld.createConstPtr();
    Lucee::FieldPtr<double> dgPtr = dgFld.createPtr();

    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
    {
      fvFld.setPtr(fvPtr, i);
      fvFld.setPtr(fvPtr_m, i-1);
      fvFld.setPtr(fvPtr_p, i+1);
      dgFld.setPtr(dgPtr, i);
// compute values at DG nodes on cell edges by simple averaging
      dgPtr[0] = 0.5*(fvPtr_m[0] + fvPtr[0]);
      dgPtr[1] = 0.5*(fvPtr[0] + fvPtr_p[0]);
    }

    if (extrapolateNodes)
    {
// modify nodes on domain boundary using extrapolation, rather than using ghost cells
      if (localRgn.getLower(0) == fvFld.getGlobalLower(0))
      { // on correct processor
        fvFld.setPtr(fvPtr, localRgn.getLower(0));
        fvFld.setPtr(fvPtr_p, localRgn.getLower(0)+1);
        dgFld.setPtr(dgPtr, localRgn.getLower(0));
// extrapolate linearly
        dgPtr[0] = 1.5*fvPtr[0] - 0.5*fvPtr_p[0];
      }
      if (localRgn.getUpper(0) == fvFld.getGlobalUpper(0))
      { // on correct processor
        fvFld.setPtr(fvPtr, localRgn.getUpper(0)-1); // getUpper() returns one beyond last index
        fvFld.setPtr(fvPtr_m, localRgn.getUpper(0)-2);
        dgFld.setPtr(dgPtr, localRgn.getUpper(0)-1);
// extrapolate linearly
        dgPtr[1] = 1.5*fvPtr[0] - 0.5*fvPtr_m[0];
      }
    }
    
    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  FiniteVolumeToLinearDGUpdater<NDIM>::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<NDIM, double>));
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class Lucee::FiniteVolumeToLinearDGUpdater<1>;
  template class Lucee::FiniteVolumeToLinearDGUpdater<2>;
  template class Lucee::FiniteVolumeToLinearDGUpdater<3>;
}
