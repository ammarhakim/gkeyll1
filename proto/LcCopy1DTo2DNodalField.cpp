/**
 * @file	LcCopy1DTo2DNodalField.cpp
 *
 * @brief	Updater to copy 1D field to 2D field.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCopy1DTo2DNodalField.h>

namespace Lucee
{
  const char *Copy1DTo2DNodalField::id = "Copy1DTo2DNodalField";

  Copy1DTo2DNodalField::Copy1DTo2DNodalField()
    : Lucee::UpdaterIfc()
  {
  }

  void
  Copy1DTo2DNodalField::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);

// get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify element to use using 'basis'");
  }

  void
  Copy1DTo2DNodalField::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  Copy1DTo2DNodalField::update(double t)
  {
    return Lucee::UpdaterStatus();
  }

  void
  Copy1DTo2DNodalField::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}
