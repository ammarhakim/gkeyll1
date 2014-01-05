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
#include <LcCenterOfMassUpdater.h>
#include <LcCopy1DTo2DNodalField.h>
#include <LcDGDiffusionUpdater1D.h>
#include <LcDistFuncMomentCalc1D.h>
#include <LcDistFuncReflectionBcUpdater.h>
#include <LcEnergyFromStreamAndVortUpdater.h>
#include <LcEnergyFromStreamFunctionUpdater.h>
#include <LcEnstrophyUpdater.h>
#include <LcIntegrateFieldAlongLine.h>
#include <LcIntegrateNodalField.h>
#include <LcLinEmGke1dHamilPertUpdater.h>
#include <LcMaxwellTm2DUpdater.h>
#include <LcModalDg1DDiffusionUpdater.h>
#include <LcModalDg1DHyperDiffusionUpdater.h>
#include <LcModalDg1DLocalDGUpdater.h>
#include <LcModalDg1DSymmetricDDGUpdater.h>
#include <LcModalDg1DUpdater.h>
#include <LcModalDgLimiter1DUpdater.h>
#include <LcModalL2NormUpdater.h>
#include <LcMusclHancock1DUpdater.h>
#include <LcNodalGradientUpdater.h>
#include <LcNodalPoissonBracketUpdater.h>
#include <LcNonLinEmGke1dHamilUpdater.h>
#include <LcNormGradPhiUpdater.h>
#include <LcProtoSolverRegistry.h>
#include <LcRecordFieldDerivInCell.h>
#include <LcRecordFieldInCell.h>
#include <LcRectSecondOrderCentralDiffUpdater.h>
#include <LcRegisteredObjList.h>
#include <LcSheathParticleSource1x1v.h>

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
      .append<Lucee::CenterOfMassUpdater<1> >()
      .append<Lucee::CenterOfMassUpdater<2> >()
      .append<Lucee::CenterOfMassUpdater<3> >()

      .append<Lucee::Copy1DTo2DNodalField>()
      .append<Lucee::DistFuncMomentCalc1D>()
      .append<Lucee::DistFuncReflectionBcUpdater>()
      .append<Lucee::EnergyFromStreamAndVortUpdater>()
      .append<Lucee::EnergyFromStreamFunctionUpdater>()
    
      .append<Lucee::EnstrophyUpdater<1> >()
      .append<Lucee::EnstrophyUpdater<2> >()
      .append<Lucee::EnstrophyUpdater<3> >()

      .append<Lucee::IntegrateFieldAlongLine<1> >()
      .append<Lucee::IntegrateFieldAlongLine<2> >()
      .append<Lucee::IntegrateFieldAlongLine<3> >()

      .append<Lucee::IntegrateNodalField<1> >()
      .append<Lucee::IntegrateNodalField<2> >()
      .append<Lucee::IntegrateNodalField<3> >()

      .append<Lucee::MaxwellTm2DUpdater>()
      .append<Lucee::ModalDg1DDiffusionUpdater>()
      .append<Lucee::ModalDg1DHyperDiffusionUpdater>()
      .append<Lucee::ModalDg1DLocalDGUpdater>()
      .append<Lucee::ModalDg1DSymmetricDDGUpdater>()
      .append<Lucee::ModalDg1DUpdater>()
      .append<Lucee::ModalDgLimiter1DUpdater>()
      .append<Lucee::ModalL2NormUpdater>()
      .append<Lucee::MusclHancock1DUpdater>()
      .append<Lucee::DGDiffusionUpdater1D>()

      .append<Lucee::NodalGradientUpdater<1> >()
      .append<Lucee::NodalGradientUpdater<2> >()
      .append<Lucee::NodalGradientUpdater<3> >()

      .append<Lucee::NodalPoissonBracketUpdater>()

      .append<Lucee::NormGradPhiUpdater<1> >()
      .append<Lucee::NormGradPhiUpdater<2> >()
      .append<Lucee::NormGradPhiUpdater<3> >()

      .append<Lucee::RecordFieldInCell<1> >()
      .append<Lucee::RecordFieldInCell<2> >()
      .append<Lucee::RecordFieldInCell<3> >()

      .append<Lucee::RecordFieldDerivInCell<1> >()
      .append<Lucee::RecordFieldDerivInCell<2> >()
      .append<Lucee::RecordFieldDerivInCell<3> >()

      .append<Lucee::RectSecondOrderCentralDiffUpdater<1> >()
      .append<Lucee::RectSecondOrderCentralDiffUpdater<2> >()
      .append<Lucee::RectSecondOrderCentralDiffUpdater<3> >()
    
      .append<Lucee::LinEmGke1dPertHamilUpdater>()
      .append<Lucee::NonLinEmGke1dHamilUpdater>()

      .append<Lucee::SheathParticleSource1x1v>();

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
