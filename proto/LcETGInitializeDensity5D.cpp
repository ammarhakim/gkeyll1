/**
 * @file	LcETGInitializeDensity5D.cpp
 *
 * @brief	Updater to compute 2d moments of a 5d distribution function with an additional weighting function.
 * Currently only works for 0th nIn
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcETGInitializeDensity5D.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *ETGInitializeDensity5D::id = "ETGInitializeDensity5D";

  ETGInitializeDensity5D::ETGInitializeDensity5D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  ETGInitializeDensity5D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("ETGInitializeDensity5D::readInput: Must specify polyOrder");

    // get hold of 5D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<5> >("basis5d"))
      nodalBasis5d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<5> >("basis5d");
    else
      throw Lucee::Except(
        "ETGInitializeDensity5D::readInput: Must specify 5D element to use using 'basis5d'");

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis3d"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis3d");
    else
      throw Lucee::Except(
        "ETGInitializeDensity5D::readInput: Must specify 3D element to use using 'basis1d'");

    // get nIn to compute
    if (tbl.hasNumber("constantDensity"))
    constantDensity = tbl.getNumber("constantDensity");
    else
      throw Lucee::Except(
        "ETGInitializeDensity5D::readInput: Must specify density using 'constantDensity'");
  }

  void
  ETGInitializeDensity5D::initialize()
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
  ETGInitializeDensity5D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<5>& grid 
      = this->getGrid<Lucee::StructuredGridBase<5> >();

    // get input field (3d)
    const Lucee::Field<3, double>& nIn = this->getInp<Lucee::Field<3, double> >(0);
    // get output field (5d)
    Lucee::Field<5, double>& distF = this->getOut<Lucee::Field<5, double> >(0);

    // local region to update (This is the 5D region. The 3D region is
    // assumed to have the same cell layout as the X-direction of the 5D region)
    Lucee::Region<5, int> localExtRgn = distF.getExtRegion();
    //Lucee::Region<5, int> localRgn = grid.getLocalRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> nPtr = nIn.createConstPtr();
    Lucee::FieldPtr<double> distFPtr = distF.createPtr();

    int idx[5];
    Lucee::RowMajorSequencer<5> seq(localExtRgn);
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    unsigned nlocal5d = nodalBasis5d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    // Loop over each cell in 5D space
    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);

      nIn.setPtr(nPtr, idx[0], idx[1], idx[2]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd scaleFactors3d(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        scaleFactors3d(i) = constantDensity/nPtr[i];

      Eigen::VectorXd scaleFactors5d = mappingMatrix*scaleFactors3d;

      for (int i = 0; i < nlocal5d; i++)
        distFPtr[i] = scaleFactors5d(i)*distFPtr[i];
    }

    return Lucee::UpdaterStatus();
  }

  void
  ETGInitializeDensity5D::declareTypes()
  {
    // computed density (3d)
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // output scaled distribution function (5d)
    this->appendOutVarType(typeid(Lucee::Field<5, double>));
  }
}
