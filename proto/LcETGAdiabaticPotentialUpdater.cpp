/**
 * @file	LcETGAdiabaticPotentialUpdater.cpp
 *
 * @brief	Updater to compute the potential for ETG test problem
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcETGAdiabaticPotentialUpdater.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *ETGAdiabaticPotentialUpdater::id = "ETGAdiabaticPotentialUpdater";

  ETGAdiabaticPotentialUpdater::ETGAdiabaticPotentialUpdater()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  ETGAdiabaticPotentialUpdater::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater::readInput: Must specify 2D element to use using 'basis'");

    if (tbl.hasNumber("kzfTimesRhoSquared"))
      kzfTimesRhoSquared = tbl.getNumber("kzfTimesRhoSquared");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater::readInput: Must specify value for 'kzfTimesRhoSquared'");

    if (tbl.hasNumber("adiabaticTemp"))
      adiabaticTemp = tbl.getNumber("adiabaticTemp");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater::readInput: Must specify value for 'adiabaticTemp' (in eV)");

    if (tbl.hasNumber("adiabaticCharge"))
      adiabaticCharge = tbl.getNumber("adiabaticCharge");
    else
      throw Lucee::Except(
        "ETGAdiabaticPotentialUpdater::readInput: Must specify value for 'adiabaticCharge'");
  }

  void
  ETGAdiabaticPotentialUpdater::initialize()
  {
    // call base class method
    Lucee::UpdaterIfc::initialize();

    // get number of nodes in 1D and 2D
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    // get volume interpolation matrices for 2d element
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    std::vector<double> volWeights2d(nVolQuad2d);
    Lucee::Matrix<double> tempVolQuad2d(nVolQuad2d, nlocal2d);
    Lucee::Matrix<double> tempVolCoords2d(nVolQuad2d, 3);

    volQuad2d.reset(nVolQuad2d, nlocal2d);

    // get data needed for Gaussian quadrature
    nodalBasis2d->getGaussQuadData(tempVolQuad2d, tempVolCoords2d, volWeights2d);
    for (int volIndex = 0; volIndex < nVolQuad2d; volIndex++)
      volQuad2d.weights(volIndex) = volWeights2d[volIndex];
    
    copyLuceeToEigen(tempVolQuad2d, volQuad2d.interpMat);
  }

  Lucee::UpdaterStatus
  ETGAdiabaticPotentialUpdater::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // get input fields (2d)
    const Lucee::Field<2, double>& nKineticIn = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& nAdiabaticIn = this->getInp<Lucee::Field<2, double> >(1);
    // get output field (2d)
    Lucee::Field<2, double>& phiOut = this->getOut<Lucee::Field<2, double> >(0);

    Lucee::Region<2, int> localRgn = grid.getLocalRegion();
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> nKineticPtr = nKineticIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nAdiabaticPtr = nAdiabaticIn.createConstPtr();
    Lucee::FieldPtr<double> phiOutPtr = phiOut.createPtr();

    int idx[2];
    Lucee::RowMajorSequencer<2> seq(localRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();

    // Find ion density in center of domain (averaged from corners of four adj. cells)
    double nAdiabaticAtCenter = 0.0;

    for (int i = 0; i < 2; i++)
      idx[i] = (globalRgn.getLower(i)+globalRgn.getUpper(i))/2;

    // NOTE: ONLY WORKS FOR POLYORDER 1
    // UPPER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx);
    nAdiabaticAtCenter += 0.25*nAdiabaticPtr[0];
    // LOWER LEFT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0], idx[1]-1);
    nAdiabaticAtCenter += 0.25*nAdiabaticPtr[2];
    // LOWER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1]-1);
    nAdiabaticAtCenter += 0.25*nAdiabaticPtr[3];
    // UPPER RIGHT
    nAdiabaticIn.setPtr(nAdiabaticPtr, idx[0]-1, idx[1]);
    nAdiabaticAtCenter += 0.25*nAdiabaticPtr[1];

    int NY_TOTAL = globalRgn.getUpper(1) - globalRgn.getLower(1);

    // Loop over each cell
    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Loop over cells in y to compute average phi
      // TODO: Change inputs in future, compute <rho>_y using a moment updater instead
      Eigen::VectorXd rhoYAvgVec = Eigen::VectorXd::Zero(nlocal2d);

      for (int iy = globalRgn.getLower(1); iy < globalRgn.getUpper(1); iy++)
      {
        idx[0] = ix;
        idx[1] = iy;

        grid.setIndex(idx);

        nKineticIn.setPtr(nKineticPtr, idx);
        nAdiabaticIn.setPtr(nAdiabaticPtr, idx);

        Eigen::VectorXd rhoVec(nlocal2d);
        for (int i = 0; i < nlocal2d; i++)
          rhoVec(i) = -(adiabaticTemp*fabs(adiabaticCharge)/(nAdiabaticAtCenter*adiabaticCharge))
            *(nKineticPtr[i] - nAdiabaticPtr[i]);

        phiOut.setPtr(phiOutPtr, idx);
        for (int i = 0; i < nlocal2d; i++)
          phiOutPtr[i] = rhoVec(i);


        // NOTE: ONLY WORKS FOR POLYORDER 1
        rhoYAvgVec(0) += (rhoVec(0) + rhoVec(2));
        rhoYAvgVec(1) += (rhoVec(1) + rhoVec(3));
        rhoYAvgVec(2) += (rhoVec(0) + rhoVec(2));
        rhoYAvgVec(3) += (rhoVec(1) + rhoVec(3));
      }

      rhoYAvgVec /= 2*NY_TOTAL;
      
      // Accumulate kzf term contribution
      if (kzfTimesRhoSquared != 1.0)
      {
        for (int iy = globalRgn.getLower(1); iy < globalRgn.getUpper(1); iy++)
        {
          idx[0] = ix;
          idx[1] = iy;
          phiOut.setPtr(phiOutPtr, idx);

          for (int i = 0; i < nlocal2d; i++)
            phiOutPtr[i] = (1 - kzfTimesRhoSquared)/kzfTimesRhoSquared*rhoYAvgVec(i) + phiOutPtr[i];
        }
      }
    }

    return Lucee::UpdaterStatus();
  }

  void
  ETGAdiabaticPotentialUpdater::declareTypes()
  {
    // inputs: kinetic density, adiabatic density
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // output potential
    this->appendOutVarType(typeid(Lucee::Field<2, double>));
  }

  void
  ETGAdiabaticPotentialUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
