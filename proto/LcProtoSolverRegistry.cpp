/**
 * @file	LcProtoSolverRegistry.cpp
 *
 * @brief	Method for registering basic solver object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCopy1DTo2DNodalField.h>
#include <LcDistFuncMomentCalc1D.h>
#include <LcEnergyFromStreamAndVortUpdater.h>
#include <LcEnergyFromStreamFunctionUpdater.h>
#include <LcEnstrophyUpdater.h>
#include <LcIntegrateNodalField.h>
#include <LcMaxwellTm2DUpdater.h>
#include <LcModalDg1DUpdater.h>
#include <LcModalDg1DDiffusionUpdater.h>
#include <LcModalDg1DHyperDiffusionUpdater.h>
#include <LcModalDg1DLocalDGUpdater.h>
#include <LcModalDg1DSymmetricDDGUpdater.h>
#include <LcModalDgLimiter1DUpdater.h>
#include <LcModalL2NormUpdater.h>
#include <LcMusclHancock1DUpdater.h>
#include <LcNodalGradientUpdater.h>
#include <LcNodalPoissonBracketUpdater.h>
#include <LcNormGradPhiUpdater.h>
#include <LcProtoSolverRegistry.h>
#include <LcRecordFieldInCell.h>
#include <LcRectSecondOrderCentralDiffUpdater.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

#ifdef HAVE_PETSC
# include <LcContFromDisContUpdater.h>
#endif

#ifdef HAVE_FFTW3
# include <LcPeriodicPoisson2DUpdater.h>
#endif

namespace Lucee
{
  void
  registerProtoSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::UpdaterIfc> >
      ::Instance()
      .append<Lucee::MaxwellTm2DUpdater>()
      .append<Lucee::MusclHancock1DUpdater>()
      .append<Lucee::RectSecondOrderCentralDiffUpdater<1> >()
      .append<Lucee::RectSecondOrderCentralDiffUpdater<2> >()
      .append<Lucee::RectSecondOrderCentralDiffUpdater<3> >()
      .append<Lucee::ModalDg1DUpdater>()
      .append<Lucee::ModalDg1DDiffusionUpdater>()
      .append<Lucee::ModalDg1DHyperDiffusionUpdater>()
      .append<Lucee::ModalDg1DLocalDGUpdater>()
      .append<Lucee::ModalDg1DSymmetricDDGUpdater>()
      .append<Lucee::ModalDgLimiter1DUpdater>()
      .append<Lucee::ModalL2NormUpdater>()
      .append<Lucee::NodalPoissonBracketUpdater>()
      .append<Lucee::EnergyFromStreamFunctionUpdater>()
      .append<Lucee::EnergyFromStreamAndVortUpdater>()
      .append<Lucee::EnstrophyUpdater>()
      .append<Lucee::NodalGradientUpdater>()
      .append<Lucee::DistFuncMomentCalc1D>()
      .append<Lucee::IntegrateNodalField<1> >()
      .append<Lucee::IntegrateNodalField<2> >()
      .append<Lucee::IntegrateNodalField<3> >()
      .append<Lucee::RecordFieldInCell<1> >()
      .append<Lucee::RecordFieldInCell<2> >()
      .append<Lucee::RecordFieldInCell<3> >()
      .append<Lucee::Copy1DTo2DNodalField>()
      .append<Lucee::NormGradPhiUpdater<1> >()
      .append<Lucee::NormGradPhiUpdater<2> >()
      .append<Lucee::NormGradPhiUpdater<3> >();
#ifdef HAVE_PETSC
      Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::UpdaterIfc> >
      ::Instance()
      .append<Lucee::ContFromDisContUpdater<1> >();
#endif
#ifdef HAVE_FFTW3
      Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::UpdaterIfc> >
      ::Instance()
      .append<Lucee::PeriodicPoisson2DUpdater>();
#endif
  }
}
