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
      throw Lucee::Except("Copy1DTo2DNodalField::readInput: Must specify element to use using 'basis'");
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
// get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

// get input field (1d)
    const Lucee::Field<1, double>& fld1d = this->getInp<Lucee::Field<1, double> >(0);
// get output field (2D)
    Lucee::Field<2, double>& fld2d = this->getOut<Lucee::Field<2, double> >(0);

    unsigned polyOrder = 1;
    if (nodalBasis->getNumNodes() == 2)
      polyOrder = 1;
    else if (nodalBasis->getNumNodes() == 3)
      polyOrder = 2;
    else
    {
      Lucee::Except lce("Copy1DTo2DNodalField::update: element with nodes ");
      lce << nodalBasis->getNumNodes() << " not supported";
      throw lce;
    }

// local region to update (This is the 2D region. The 1D region is
// assumed to have the same cell layout as the X-direction of the 2D region)
    Lucee::Region<2, int> localRgn = grid.getLocalRegion();

    Lucee::FieldPtr<double> fld2dPtr = fld2d.createPtr();
    std::vector<double> data1d(nodalBasis->getNumNodes());
// loop over all X-direction cells
    for (int i=localRgn.getLower(0)-1; i<localRgn.getUpper(0); ++i)
    {
      nodalBasis->setIndex(i);
// extract data from current cell
      nodalBasis->extractFromField(fld1d, data1d);

// copy this into alll Y-direction cells
      for (int j=localRgn.getLower(1); j<localRgn.getUpper(1)+1; ++j)
      {
        fld2d.setPtr(fld2dPtr, i, j);
// copy data based on polynomial order
        if (polyOrder == 1)
          fld2dPtr[0] = data1d[0];
        else
        {
          fld2dPtr[0] = data1d[0];
          fld2dPtr[1] = data1d[1];
          fld2dPtr[2] = data1d[0];
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
