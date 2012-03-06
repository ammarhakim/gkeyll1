/**
 * @file	LcFemPoisson1DUpdater.cpp
 *
 * @brief	Updater to solve Poisson equations with FEM scheme in 1D.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFemPoisson1DUpdater.h>
#include <LcStructGridField.h>

namespace Lucee
{
  const char *FemPoisson1DUpdater::id = "FemPoisson1D";

  FemPoisson1DUpdater::FemPoisson1DUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  FemPoisson1DUpdater::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc>("nodalEement"))
      nodalFem = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc>("nodalElement");
    else
      throw Lucee::Except("FemPoisson1DUpdater::readInput: Must specify element to use using 'nodalElement'");
  }

  void
  FemPoisson1DUpdater::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  FemPoisson1DUpdater::update(double t)
  {

    return Lucee::UpdaterStatus();
  }

  void
  FemPoisson1DUpdater::declareTypes()
  {
// takes one input (source terms)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
// returns one output, solution
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
