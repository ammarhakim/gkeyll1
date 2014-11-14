/**
 * @file	LcETGAdiabaticPotentialUpdater3D.cpp
 *
 * @brief	Updater to compute the potential for ETG test problem
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcETGAdiabaticPotentialUpdater3D.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *ETGAdiabaticPotentialUpdater3D::id = "ETGAdiabaticPotentialUpdater3D";

  ETGAdiabaticPotentialUpdater3D::ETGAdiabaticPotentialUpdater3D()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  ETGAdiabaticPotentialUpdater3D::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 3D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis3d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater3D::readInput: Must specify 3D element to use using 'basis'");

    if (tbl.hasNumber("adiabaticTemp"))
      adiabaticTemp = tbl.getNumber("adiabaticTemp");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater3D::readInput: Must specify value for 'adiabaticTemp' (in eV)");

    if (tbl.hasNumber("adiabaticCharge"))
      adiabaticCharge = tbl.getNumber("adiabaticCharge");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater3D::readInput: Must specify value for 'adiabaticCharge'");
  }

  void
  ETGAdiabaticPotentialUpdater3D::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get number of nodes in 3D
    unsigned nlocal3d = nodalBasis3d->getNumNodes();

    // get volume interpolation matrices for 3d element
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();
    std::vector<double> volWeights3d(nVolQuad3d);
    Lucee::Matrix<double> tempVolQuad3d(nVolQuad3d, nlocal3d);
    Lucee::Matrix<double> tempVolCoords3d(nVolQuad3d, 3);

    volQuad3d.reset(nVolQuad3d, nlocal3d);

    // get data needed for Gaussian quadrature
    nodalBasis3d->getGaussQuadData(tempVolQuad3d, tempVolCoords3d, volWeights3d);
    for (int volIndex = 0; volIndex < nVolQuad3d; volIndex++)
      volQuad3d.weights(volIndex) = volWeights3d[volIndex];
    
    copyLuceeToEigen(tempVolQuad3d, volQuad3d.interpMat);
  }

  Lucee::UpdaterStatus
  ETGAdiabaticPotentialUpdater3D::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // get input fields (3d)
    const Lucee::Field<3, double>& nKineticIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& nAdiabaticIn = this->getInp<Lucee::Field<3, double> >(1);
    // get output field (3d)
    Lucee::Field<3, double>& phiOut = this->getOut<Lucee::Field<3, double> >(0);

    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::Region<3, int> globalRgn = grid.getGlobalRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> nKineticPtr = nKineticIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nAdiabaticPtr = nAdiabaticIn.createConstPtr();
    Lucee::FieldPtr<double> phiOutPtr = phiOut.createPtr();

    int idx[3];
    Lucee::RowMajorSequencer<3> seq(localRgn);
    unsigned nlocal3d = nodalBasis3d->getNumNodes();
    int nVolQuad3d = nodalBasis3d->getNumGaussNodes();

    // Find ion density in center of domain (averaged from corners of four adj. cells)
    double nAdiabaticAtCenter = 0.0;

    for (int i = 0; i < 3; i++)
      idx[i] = (globalRgn.getLower(i)+globalRgn.getUpper(i))/2;

    // NOTE: ONLY WORKS FOR POLYORDER 1
    // UPPER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[0];
    // LOWER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0], idx[1]-1, idx[2]);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[2];
    // LOWER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1]-1, idx[2]);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[3];
    // UPPER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1], idx[2]);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[1];
    // UPPER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0], idx[1], idx[2]-1);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[4];
    // LOWER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0], idx[1]-1, idx[2]-1);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[6];
    // LOWER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1]-1, idx[2]-1);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[7];
    // UPPER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1], idx[2]-1);
    nAdiabaticAtCenter += 1/8.0*nAdiabaticPtr[5];

    while (seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);

      nKineticIn.setPtr(nKineticPtr, idx);
      nAdiabaticIn.setPtr(nAdiabaticPtr, idx);

      Eigen::VectorXd rhoVec(nlocal3d);
      for (int i = 0; i < nlocal3d; i++)
        rhoVec(i) = -(adiabaticTemp*fabs(adiabaticCharge)/(nAdiabaticAtCenter*adiabaticCharge))
          *(nKineticPtr[i] - nAdiabaticPtr[i]);

      phiOut.setPtr(phiOutPtr, idx);
      for (int i = 0; i < nlocal3d; i++)
        phiOutPtr[i] = rhoVec(i);
    }

    return Lucee::UpdaterStatus();
  }

  void
  ETGAdiabaticPotentialUpdater3D::declareTypes()
  {
    // inputs: kinetic density, adiabatic density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // output potential
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }

  void
  ETGAdiabaticPotentialUpdater3D::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
