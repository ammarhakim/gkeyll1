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
// NOTE: This method only works for polyOrder 1 or 2 serendipity
// elements. Eventually, I need to fix this, but not sure how. (Ammar
// Hakim, August 2012).

// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get input field (1d)
    const Lucee::Field<1, double>& fld1d = this->getInp<Lucee::Field<1, double> >(0);
// get output field (2D)
    Lucee::Field<2, double>& fld2d = this->getOut<Lucee::Field<2, double> >(0);

    unsigned polyOrder = 1;
    if (fld1d.getNumComponents() == 1)
      polyOrder = 1;
    else if (fld1d.getNumComponents() == 2)
      polyOrder = 2;
    else
    {
      Lucee::Except lce("Copy1DTo2DNodalField::update: element not supported");
      throw lce;
    }

// local region to update (This is the 2D region. The 1D region is
// assumed to have the same cell layout as the X-direction of the 2D region)
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> fld1dPtr = fld1d.createConstPtr();
    Lucee::FieldPtr<double> fld2dPtr = fld2d.createPtr();

// loop over all X-direction cells
    for (int i=localRgn.getLower(0)-1; i<localRgn.getUpper(0)+1; ++i)
    {
      fld1d.setPtr(fld1dPtr, i);
// copy this into all Y-direction cells
      for (int j=localRgn.getLower(1)-1; j<localRgn.getUpper(1)+1; ++j)
      {
        fld2d.setPtr(fld2dPtr, i, j);
// copy data based on polynomial order
        if (polyOrder == 1)
          fld2dPtr[0] = fld1dPtr[0];
        else
        {
          fld2dPtr[0] = fld1dPtr[0];
          fld2dPtr[1] = fld1dPtr[1];
          fld2dPtr[2] = fld1dPtr[0];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  Copy1DTo2DNodalField::declareTypes()
  {
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }
}
