/**
 * @file	LcDGExplicitIncrementFiveMomentSrcUpdater.cpp
 *
 * @brief	Calculates the increment for explicitly updating 5-moment source terms using a discontinuous Galerkin method
 */

// lucee includes
#include <LcDGExplicitIncrementFiveMomentSrcUpdater.h>
#include <LcStructGridField.h>

// eigen inlcudes
#include <Eigen/Eigen>

namespace Lucee
{

  static const unsigned X = 0;
  static const unsigned Y = 1;
  static const unsigned Z = 2;

  static const unsigned RHO = 0;
  static const unsigned RHOUX = 1;
  static const unsigned RHOUY = 2;
  static const unsigned RHOUZ = 3;
  static const unsigned ER = 4;

  static const unsigned EX = 0;
  static const unsigned EY = 1;
  static const unsigned EZ = 2;
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;
  static const unsigned PHIE = 6;
  static const unsigned PHIM = 7;

// number of components for EM (8) and fluid (5) variables
  static const unsigned numEMcomp = 8;
  static const unsigned numFldcomp = 5;

// set ids for module system
  template <> const char *DGExplicitIncrementFiveMomentSrcUpdater<1>::id = "DGExplicitIncrementFiveMomentSrc1D";
  template <> const char *DGExplicitIncrementFiveMomentSrcUpdater<2>::id = "DGExplicitIncrementFiveMomentSrc2D";
  template <> const char *DGExplicitIncrementFiveMomentSrcUpdater<3>::id = "DGExplicitIncrementFiveMomentSrc3D";

  template <unsigned NDIM>
  void
  DGExplicitIncrementFiveMomentSrcUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      basis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("DGExplicitIncrementFiveMomentSrcUpdater::readInput: Must specify configuration-space basis using 'basis'");

    nFluids = (int) tbl.getNumber("numFluids");

// charge on each fluid species
    charge = tbl.getNumVec("charge");
    if (charge.size() != nFluids)
    {
      Lucee::Except lce("DGExplicitIncrementFiveMomentSrcUpdater::readInput: 'charge' table should have exactly ");
      lce << nFluids << " entries. Provided " << charge.size() << " instead.";
      throw lce;
    }

// mass of each fluid species
    mass = tbl.getNumVec("mass");
    if (mass.size() != nFluids)
    {
      Lucee::Except lce("DGExplicitIncrementFiveMomentSrcUpdater::readInput: 'mass' table should have exactly ");
      lce << nFluids << " entries. Provided " << mass.size() << " instead.";
      throw lce;
    }

// permittivity of free space
    epsilon0 = tbl.getNumber("epsilon0");

    qbym.resize(nFluids);
    qbym2.resize(nFluids);
    for (unsigned i=0; i<nFluids; ++i)
    {
      qbym[i] = charge[i]/mass[i];
      qbym2[i] = qbym[i]*qbym[i];
    }
  }

  template <unsigned NDIM>
  void
  DGExplicitIncrementFiveMomentSrcUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  DGExplicitIncrementFiveMomentSrcUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();

    unsigned nlocal = basis->getNumNodes();

    double dt = t-this->getCurrTime();

// get pointers to each input fluid
    std::vector<const Lucee::Field<NDIM, double>* > fluidsIn;
    for (unsigned i=0; i<nFluids; ++i)
      fluidsIn.push_back(&this->getInp<Lucee::Field<NDIM, double> >(i));

// get input EM field (this is last out parameter)
    const Lucee::Field<NDIM, double>& emFieldIn = this->getInp<Lucee::Field<NDIM, double> >(nFluids);
    
// get pointers to each output fluid
    std::vector<Lucee::Field<NDIM, double>* > fluidsNew;
    for (unsigned i=0; i<nFluids; ++i)
      fluidsNew.push_back(&this->getOut<Lucee::Field<NDIM, double> >(i));

// get output EM field (this is last out parameter)
    Lucee::Field<NDIM, double>& emFieldNew = this->getOut<Lucee::Field<NDIM, double> >(nFluids);

    Lucee::ConstFieldPtr<double> fPtrIn = fluidsIn[0]->createConstPtr();
    Lucee::ConstFieldPtr<double> emPtrIn = emFieldIn.createConstPtr();
    Lucee::FieldPtr<double> fPtrNew = fluidsNew[0]->createPtr();
    Lucee::FieldPtr<double> emPtrNew = emFieldNew.createPtr();

    int idx[NDIM];
    Lucee::Region<NDIM, int> localRgn = emFieldNew.getRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);

      emFieldIn.setPtr(emPtrIn, idx);
      emFieldNew.setPtr(emPtrNew, idx);
      Eigen::Vector3d E;
      for (unsigned i=0; i<nlocal; ++i)
      {

        double keold = 0.0;
        E(0) = 0.0; E(1) = 0.0; E(2) = 0.0;
// update solution for fluids (solution is at half-time step)
        for (unsigned n=0; n<nFluids; ++n)
        {
          fluidsIn[n]->setPtr(fPtrIn, idx);
          fluidsNew[n]->setPtr(fPtrNew, idx);

          double wp2 = fPtrIn[i*numFldcomp+RHO]*qbym2[n]/epsilon0;
          E(0) += -qbym[n]/epsilon0*fPtrIn[i*numFldcomp+RHOUX];
          E(1) += -qbym[n]/epsilon0*fPtrIn[i*numFldcomp+RHOUY]; 
          E(2) += -qbym[n]/epsilon0*fPtrIn[i*numFldcomp+RHOUZ];

// compute old kinetic energy before momenta are over-written
          keold = 0.5*(fPtrIn[i*numFldcomp+RHOUX]*fPtrIn[i*numFldcomp+RHOUX] + fPtrIn[i*numFldcomp+RHOUY]*fPtrIn[i*numFldcomp+RHOUY] + fPtrIn[i*numFldcomp+RHOUZ]*fPtrIn[i*numFldcomp+RHOUZ])/fPtrIn[i*numFldcomp+RHO];

          fPtrNew[i*numFldcomp+RHO] = 0.0;
          fPtrNew[i*numFldcomp+RHOUX] = (wp2*epsilon0*emPtrIn[i*numEMcomp+EX] 
             + qbym2[n]*(fPtrIn[i*numFldcomp+RHOUY]*emPtrIn[i*numEMcomp+BZ] - fPtrIn[i*numFldcomp+RHOUZ]*emPtrIn[i*numEMcomp+BY]))/qbym[n];
          fPtrNew[i*numFldcomp+RHOUY] = (wp2*epsilon0*emPtrIn[i*numEMcomp+EY] 
             + qbym2[n]*(fPtrIn[i*numFldcomp+RHOUZ]*emPtrIn[i*numEMcomp+BX] - fPtrIn[i*numFldcomp+RHOUX]*emPtrIn[i*numEMcomp+BZ]))/qbym[n];
          fPtrNew[i*numFldcomp+RHOUZ] = (wp2*epsilon0*emPtrIn[i*numEMcomp+EZ] 
             + qbym2[n]*(fPtrIn[i*numFldcomp+RHOUX]*emPtrIn[i*numEMcomp+BY] - fPtrIn[i*numFldcomp+RHOUY]*emPtrIn[i*numEMcomp+BX]))/qbym[n];
// energy equation (there is no explicit energy source, so just
// recompute new kinetic energy to update total energy)
          fPtrNew[i*numFldcomp+ER] = fPtrIn[i*numFldcomp+ER] - keold 
             + 0.5*(fPtrNew[i*numFldcomp+RHOUX]*fPtrNew[i*numFldcomp+RHOUX] 
             + fPtrNew[i*numFldcomp+RHOUY]*fPtrNew[i*numFldcomp+RHOUY] 
             + fPtrNew[i*numFldcomp+RHOUZ]*fPtrNew[i*numFldcomp+RHOUZ])/fPtrIn[i*numFldcomp+RHO];
        }
        emPtrNew[i*numEMcomp+EX] = E(0);
        emPtrNew[i*numEMcomp+EY] = E(1);
        emPtrNew[i*numEMcomp+EZ] = E(2);
      }
    }
    
    return Lucee::UpdaterStatus();
  }
  

  template <unsigned NDIM>
  void
  DGExplicitIncrementFiveMomentSrcUpdater<NDIM>::declareTypes()
  {
    this->setLastOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class DGExplicitIncrementFiveMomentSrcUpdater<1>;
  template class DGExplicitIncrementFiveMomentSrcUpdater<2>;
  template class DGExplicitIncrementFiveMomentSrcUpdater<3>;
}
