/**
 * @file	LcLinEmGke1dPertHamilUpdater.cpp
 *
 * @brief	Compute linearized perturbed Hamiltonian
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLinEmGke1dHamilPertUpdater.h>

namespace Lucee
{
  const char *LinEmGke1dPertHamilUpdater::id = "LinEmGke1DPertHamil";

  LinEmGke1dPertHamilUpdater::LinEmGke1dPertHamilUpdater()
    : UpdaterIfc()
  {
  }

  void
  LinEmGke1dPertHamilUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);
    charge = tbl.getNumber("charge");
    mass = tbl.getNumber("mass");
  }

  void
  LinEmGke1dPertHamilUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  LinEmGke1dPertHamilUpdater::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
// get input/output fields
    const Lucee::Field<2, double>& phi = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& Apar = this->getInp<Lucee::Field<2, double> >(1);
    Lucee::Field<2, double>& hamil = this->getOut<Lucee::Field<2, double> >(0);

    return Lucee::UpdaterStatus();
  }

  void
  LinEmGke1dPertHamilUpdater::declareTypes()
  {
// takes two inputs (phi, Apar)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
// returns one output, perturbed Hamiltonian
    this->appendOutVarType(typeid(Lucee::Field<2, double>));    
  }
}
