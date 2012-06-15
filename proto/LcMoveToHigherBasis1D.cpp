/**
 * @file	LcMoveToHigherBasis1D.cpp
 *
 * @brief	Updater to copy to a higher-order basis.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMoveToHigherBasis1D.h>

namespace Lucee
{
  const char *MoveToHigherBasis1D::id = "MoveToHigherBasis1D";

  MoveToHigherBasis1D::MoveToHigherBasis1D()
    : Lucee::UpdaterIfc()
  {
  }

  void
  MoveToHigherBasis1D::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::UpdaterIfc::readInput(tbl);
  }

  void
  MoveToHigherBasis1D::initialize()
  {
// call base class method
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  MoveToHigherBasis1D::update(double t)
  {
// get hold of grid
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

// get input field 
    const Lucee::Field<1, double>& low = this->getInp<Lucee::Field<1, double> >(0);
// get output field
    Lucee::Field<1, double>& high = this->getOut<Lucee::Field<1, double> >(0);

// pointers
    Lucee::ConstFieldPtr<double> lowPtr = low.createConstPtr();
    Lucee::FieldPtr<double> highPtr = high.createPtr();

// local region to update
    Lucee::Region<1, int> localRgn = low.getExtRegion();

    for (int i=localRgn.getLower(0); i<localRgn.getUpper(0); ++i)
    {
      low.setPtr(lowPtr, i);
      high.setPtr(highPtr, i);

      highPtr[0] = lowPtr[0];
      highPtr[1] = 0.5*(lowPtr[0] + lowPtr[1]);
      highPtr[2] = lowPtr[1];
    }

    return Lucee::UpdaterStatus();
  }

  void
  MoveToHigherBasis1D::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<1, double>));
  }
}
