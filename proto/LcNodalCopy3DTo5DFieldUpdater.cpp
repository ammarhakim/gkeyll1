/**
 * @file	LcNodalCopy3DTo5DFieldUpdater.cpp
 *
 * @brief	Updater to copy data from a 3D field onto a 5D field using the same basis functions.
 * Currently only works assuming 3D field lies on the (0,1,2) surface
 * Currently hard-coded to assume Serendipity element
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNodalCopy3DTo5DFieldUpdater.h>

namespace Lucee
{
  const char *NodalCopy3DTo5DFieldUpdater::id = "NodalCopy3DTo5DFieldUpdater";

  NodalCopy3DTo5DFieldUpdater::NodalCopy3DTo5DFieldUpdater()
    : Lucee::UpdaterIfc()
  {
  }

  void
  NodalCopy3DTo5DFieldUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("NodalCopy3DTo5DFieldUpdater::readInput: Must specify polyOrder");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except("NodalCopy3DTo5DFieldUpdater::readInput: Must specify element to use using 'basis3d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except("NodalCopy3DTo5DFieldUpdater::readInput: Must specify element to use using 'basis5d'");
  }

  void
  NodalCopy3DTo5DFieldUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    if (polyOrder == 1)
    {
      mappingMatrix = Eigen::MatrixXd(32,8);
      mappingMatrix << 1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1,
                       1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1,
                       1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1,
                       1,0,0,0,0,0,0,0,
                       0,1,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,1;
    }
  }

  Lucee::UpdaterStatus
  NodalCopy3DTo5DFieldUpdater::update(double t)
  {
    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get input field (3d)
    const Lucee::Field<3, double>& fld3d = this->getInp<Lucee::Field<3, double> >(0);
    // get output field (5d)
    Lucee::Field<5, double>& fld5d = this->getOut<Lucee::Field<5, double> >(0);

    // local exterior region to update (This is the 5D region)
    Lucee::Region<5, int> localExtRgn = fld5d.getExtRegion();

    Lucee::ConstFieldPtr<double> fld3dPtr = fld3d.createConstPtr();
    Lucee::FieldPtr<double> fld5dPtr = fld5d.createPtr();

    Lucee::RowMajorSequencer<5> seq(localExtRgn);
    int idx[5];

    int nlocal3d = nodalBasis3d->getNumNodes();
    int nlocal5d = nodalBasis5d->getNumNodes();

    while(seq.step())
    {
      seq.fillWithIndex(idx);
      fld3d.setPtr(fld3dPtr, idx[0], idx[1], idx[2]);
      fld5d.setPtr(fld5dPtr, idx);

      Eigen::VectorXd lowerDimVec(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        lowerDimVec(i) = fld3dPtr[i];

      // Use mapping to copy 3d data to 5d
      Eigen::VectorXd higherDimVec = mappingMatrix*lowerDimVec;

      // Store 4d data into output field
      for (int i = 0; i < nlocal5d; i++)
        fld5dPtr[i] = higherDimVec(i);
    }
 
    return Lucee::UpdaterStatus();
  }

  void
  NodalCopy3DTo5DFieldUpdater::declareTypes()
  {
    // Input: Lower-dimension field
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Higher-dimension field
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }
}
