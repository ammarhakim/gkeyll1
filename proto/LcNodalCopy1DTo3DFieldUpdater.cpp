/**
 * @file	LcNodalCopy1DTo3DFieldUpdater.cpp
 *
 * @brief	Updater to copy data from a 1D field onto a 3D field using the same basis functions.
 * Currently only works assuming 1D field lies along the 0 direction
 * Currently hard-coded to assume Serendipity element
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalCopy1DTo3DFieldUpdater.h>

namespace Lucee
{
  const char *NodalCopy1DTo3DFieldUpdater::id = "NodalCopy1DTo3DFieldUpdater";

  NodalCopy1DTo3DFieldUpdater::NodalCopy1DTo3DFieldUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  NodalCopy1DTo3DFieldUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("NodalCopy1DTo3DFieldUpdater::readInput: Must specify polyOrder");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis1d"))
      nodalBasis1d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis1d");
    else
      throw Lucee::Except("NodalCopy1DTo3DFieldUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("NodalCopy1DTo3DFieldUpdater::readInput: Must specify element to use using 'basis3d'");
  }

  void
  NodalCopy1DTo3DFieldUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    if (polyOrder == 1)
    {
      mappingMatrix = Eigen::MatrixXd(8,2);
      mappingMatrix << 1,0,
                       0,1,
                       1,0,
                       0,1,
                       1,0,
                       0,1,
                       1,0,
                       0,1;
    }
    else if (polyOrder == 2)
    {
      mappingMatrix = Eigen::MatrixXd(20,3);
      mappingMatrix << 1,0,0,
                       0,1,0,
                       0,0,1,
                       1,0,0,
                       0,0,1,
                       1,0,0,
                       0,1,0,
                       0,0,1,
                       1,0,0,
                       0,0,1,
                       1,0,0,
                       0,0,1,
                       1,0,0,
                       0,1,0,
                       0,0,1,
                       1,0,0,
                       0,0,1,
                       1,0,0,
                       0,1,0,
                       0,0,1;
    }
  }

  Lucee::UpdaterStatus
  NodalCopy1DTo3DFieldUpdater::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // get input field (1d)
    const Lucee::Field<1, double>& fld1d = this->getInp<Lucee::Field<1, double> >(0);
    // get output field (3d)
    Lucee::Field<3, double>& fld3d = this->getOut<Lucee::Field<3, double> >(0);

    // local exterior region to update (This is the 3d region)
    Lucee::Region<3, int> localExtRgn = fld3d.getExtRegion();

    Lucee::ConstFieldPtr<double> fld1dPtr = fld1d.createConstPtr();
    Lucee::FieldPtr<double> fld3dPtr = fld3d.createPtr();

    Lucee::RowMajorSequencer<3> seq(localExtRgn);
    int idx[3];

    int nlocal1d = nodalBasis1d->getNumNodes();
    int nlocal3d = nodalBasis3d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      fld1d.setPtr(fld1dPtr, idx[0]);
      fld3d.setPtr(fld3dPtr, idx);

      Eigen::VectorXd lowerDimVec(nlocal1d);
      for (int i = 0; i < nlocal1d; i++)
        lowerDimVec(i) = fld1dPtr[i];

      // Use mapping to copy 1d data to 3d
      Eigen::VectorXd higherDimVec = mappingMatrix*lowerDimVec;

      // Store 4d data into output field
      for (int i = 0; i < nlocal3d; i++)
        fld3dPtr[i] = higherDimVec(i);
    }
 
    return Lucee::UpdaterStatus();
  }

  void
  NodalCopy1DTo3DFieldUpdater::declareTypes()
  {
    // Input: Lower-dimension field
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // Output: Higher-dimension field
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
