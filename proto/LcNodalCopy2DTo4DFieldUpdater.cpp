/**
 * @file	LcNodalCopy2DTo4DFieldUpdater.cpp
 *
 * @brief	Updater to copy data from a 2D field onto a 4D field using the same basis functions.
 * Currently only works assuming 2D field lies on the (0,1) surface
 * Currently hard-coded to assume Serendipity element
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalCopy2DTo4DFieldUpdater.h>

namespace Lucee
{
  const char *NodalCopy2DTo4DFieldUpdater::id = "NodalCopy2DTo4DFieldUpdater";

  NodalCopy2DTo4DFieldUpdater::NodalCopy2DTo4DFieldUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  NodalCopy2DTo4DFieldUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("NodalCopy2DTo4DFieldUpdater::readInput: Must specify polyOrder");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("NodalCopy2DTo4DFieldUpdater::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except("NodalCopy2DTo4DFieldUpdater::readInput: Must specify element to use using 'basis4d'");
  }

  void
  NodalCopy2DTo4DFieldUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    if (polyOrder == 1)
      mappingMatrix << 1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,1,
                       1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,1,
                       1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,1,
                       1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,1;
    else if (polyOrder == 2)
      mappingMatrix << 1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1,
                       1,0,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,0,1,0,
                       1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1;
  }

  Lucee::UpdaterStatus
  NodalCopy2DTo4DFieldUpdater::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // get input field (2d)
    const Lucee::Field<2, double>& fld2d = this->getInp<Lucee::Field<2, double> >(0);
    // get output field (4D)
    Lucee::Field<4, double>& fld4d = this->getOut<Lucee::Field<4, double> >(0);

    // local region to update (This is the 4D region)
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    Lucee::ConstFieldPtr<double> fld2dPtr = fld2d.createConstPtr();
    Lucee::FieldPtr<double> fld4dPtr = fld4d.createPtr();

    Lucee::RowMajorSequencer<4> seq(localRgn);
    int idx[4];

    int nlocal2d = nodalBasis2d->getNumNodes();
    int nlocal4d = nodalBasis4d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      fld2d.setPtr(fld2dPtr, idx[0], idx[1]);
      fld4d.setPtr(fld4dPtr, idx);

      Eigen::VectorXd lowerDimVec(nlocal2d);
      for (int i = 0; i < nlocal2d; i++)
        lowerDimVec(i) = fld2dPtr[i];

      // Use mapping to copy 2d data to 4d
      Eigen::VectorXd higherDimVec = mappingMatrix*lowerDimVec;

      // Store 4d data into output field
      for (int i = 0; i < nlocal4d; i++)
        fld4dPtr[i] = higherDimVec(i);
    }
 
    return Lucee::UpdaterStatus();
  }

  void
  NodalCopy2DTo4DFieldUpdater::declareTypes()
  {
    // Input: Lower-dimension field
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Output: Higher-dimension field
    this->appendOutVarType(typeid(Lucee::Field<4, double>));
  }
}
