/**
 * @file	LcKineticTotalEnergyUpdater.cpp
 *
 * @brief	Updater to take in various DynVectors and combine them to give energy quantities.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcKineticTotalEnergyUpdater.h>
// for cutoff velocities
#include <LcDynVector.h>

namespace Lucee
{
// set id for module system
  const char *KineticTotalEnergyUpdater::id = "KineticTotalEnergyUpdater";

  KineticTotalEnergyUpdater::KineticTotalEnergyUpdater()
    : UpdaterIfc()
  {
  }

  void 
  KineticTotalEnergyUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("electronMass"))
      electronMass = tbl.getNumber("electronMass");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify electronMass");
  }

  void 
  KineticTotalEnergyUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  KineticTotalEnergyUpdater::update(double t)
  {
    // Total heat fluxes at left and right edges
    const Lucee::DynVector<double>& heatFluxesIn = this->getInp<Lucee::DynVector<double> >(0);
    // Energies computed from discrete hamiltonian
    const Lucee::DynVector<double>& energyElcIn = this->getInp<Lucee::DynVector<double> >(1);
    const Lucee::DynVector<double>& energyIonIn = this->getInp<Lucee::DynVector<double> >(2);
    // Source powers computed from discrete hamiltonian
    const Lucee::DynVector<double>& powerSrcElcIn = this->getInp<Lucee::DynVector<double> >(3);
    const Lucee::DynVector<double>& powerSrcIonIn = this->getInp<Lucee::DynVector<double> >(4);
    // Output
    Lucee::DynVector<double>& totalEnergyOut = this->getOut<Lucee::DynVector<double> >(0);

    std::vector<double> heatFluxes = heatFluxesIn.getLastInsertedData();
    std::vector<double> energyElc = energyElcIn.getLastInsertedData();
    std::vector<double> energyIon = energyIonIn.getLastInsertedData();
    std::vector<double> powerSrcElc = powerSrcElcIn.getLastInsertedData();
    std::vector<double> powerSrcIon = powerSrcIonIn.getLastInsertedData();
    
    std::vector<double> data(3);
    // Total energy at this time step
    data[0] = energyIon[0] + energyElc[0];
    // Power lost at this time step
    data[1] = heatFluxes[0] + heatFluxes[3];
    // Power added at this time step
    data[2] = powerSrcIon[0] + powerSrcElc[0];

    // Put data into the DynVector
    totalEnergyOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  KineticTotalEnergyUpdater::declareTypes()
  {
    // heat fluxes at edges
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // hamiltonian-computed energies of electrons and ions
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // hamiltonian-computed energies of SOURCE electrons and ions
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: dynvector containing energy at time t
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
}
