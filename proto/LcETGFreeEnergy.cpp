/**
 * @file	LcETGFreeEnergy.cpp
 *
 * @brief	Updater to compute the free energy for ETG problem.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathPhysConstants.h>
#include <LcETGFreeEnergy.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *ETGFreeEnergy::id = "ETGFreeEnergy";

  ETGFreeEnergy::ETGFreeEnergy()
    : Lucee::UpdaterIfc()
  {
  }

  void
  ETGFreeEnergy::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    // get hold of element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except("ETGFreeEnergy::readInput: Must specify element to use using 'basis4d'");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except("ETGFreeEnergy::readInput: Must specify element to use using 'basis2d'");

    if (tbl.hasNumber("adiabaticTemp"))
      bgAdiabaticTemp = tbl.getNumber("adiabaticTemp");
    else
      throw Lucee::Except(
        "ETGFreeEnergy::readInput: Must specify temperature of adiabatic species using 'adiabaticTemp'");

    if (tbl.hasNumber("kineticMass"))
      kineticMass = tbl.getNumber("kineticMass");
    else
      throw Lucee::Except(
        "ETGFreeEnergy::readInput: Must specify mass of kinetic species using 'kineticMass'");
  }

  void
  ETGFreeEnergy::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // local region to update
    Lucee::Region<4, int> localRgn = grid.getLocalRegion();

    // set index to first location in grid (this is okay as in this
    // updater we are assuming grid is uniform)
    Lucee::RowMajorSequencer<4> seq(localRgn);
    seq.step();
    int idx[4];
    seq.fillWithIndex(idx);
    nodalBasis4d->setIndex(idx);
    nodalBasis2d->setIndex(idx[0], idx[1]);
    
    int nlocal2d = nodalBasis2d->getNumNodes();
    int nlocal4d = nodalBasis4d->getNumNodes();

    // get data needed for Gaussian quadrature (4D)
    int nVolQuad4d = nodalBasis4d->getNumGaussNodes();
    std::vector<double> volWeights4d(nVolQuad4d);
    Lucee::Matrix<double> tempVolQuad4d(nVolQuad4d, nlocal4d);
    Lucee::Matrix<double> tempVolCoords4d(nVolQuad4d, 4);
    volQuad4d.reset(nVolQuad4d, nlocal4d);

    nodalBasis4d->getGaussQuadData(tempVolQuad4d, tempVolCoords4d, volWeights4d);
    for (int volIndex = 0; volIndex < nVolQuad4d; volIndex++)
      volQuad4d.weights(volIndex) = volWeights4d[volIndex];
    
    copyLuceeToEigen(tempVolQuad4d, volQuad4d.interpMat);

    // get data needed for Gaussian quadrature (2D)
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    std::vector<double> volWeights2d(nVolQuad2d);
    Lucee::Matrix<double> tempVolQuad2d(nVolQuad2d, nlocal2d);
    Lucee::Matrix<double> tempVolCoords2d(nVolQuad2d, 3);
    volQuad2d.reset(nVolQuad2d, nlocal2d);

    nodalBasis2d->getGaussQuadData(tempVolQuad2d, tempVolCoords2d, volWeights2d);
    for (int volIndex = 0; volIndex < nVolQuad2d; volIndex++)
      volQuad2d.weights(volIndex) = volWeights2d[volIndex];
    
    copyLuceeToEigen(tempVolQuad2d, volQuad2d.interpMat);
  }

  Lucee::UpdaterStatus
  ETGFreeEnergy::update(double t)
  {
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();
    const Lucee::Field<2, double>& bFieldIn = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& bgKineticTempIn = this->getInp<Lucee::Field<2, double> >(1);
    const Lucee::Field<2, double>& bgKineticNumDensityIn = this->getInp<Lucee::Field<2, double> >(2);
    const Lucee::Field<2, double>& perturbedKineticNumDensityIn = this->getInp<Lucee::Field<2, double> >(3);
    const Lucee::Field<4, double>& bgDistFIn = this->getInp<Lucee::Field<4, double> >(4);
    const Lucee::Field<4, double>& perturbedFIn  = this->getInp<Lucee::Field<4, double> >(5);

    Lucee::DynVector<double>& fldInt = this->getOut<Lucee::DynVector<double> >(0);

    // iterators into fields
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bgKineticTempPtr = bgKineticTempIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bgKineticNumDensityPtr = bgKineticNumDensityIn.createConstPtr();
    Lucee::ConstFieldPtr<double> perturbedKineticNumDensityPtr = perturbedKineticNumDensityIn.createConstPtr();
    Lucee::ConstFieldPtr<double> bgDistFPtr = bgDistFIn.createConstPtr();
    Lucee::ConstFieldPtr<double> perturbedFPtr = perturbedFIn.createConstPtr();

    Lucee::Region<4, int> localRgn = grid.getLocalRegion();
    Lucee::Region<4, int> globalRgn = grid.getGlobalRegion();
    Lucee::RowMajorSequencer<4> seq(localRgn);

    unsigned nlocal4d = nodalBasis4d->getNumNodes();
    unsigned nlocal2d = nodalBasis2d->getNumNodes();

    int nVolQuad4d = nodalBasis4d->getNumGaussNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    int idx[4];

    double LX = (globalRgn.getUpper(0) - globalRgn.getLower(0))*grid.getDx(0);
    double LY = (globalRgn.getUpper(1) - globalRgn.getLower(1))*grid.getDx(1);

    double localInt = 0.0;
    // loop, performing integration
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      // set index into element basis
      nodalBasis4d->setIndex(idx);
      nodalBasis2d->setIndex(idx[0], idx[1]);

      bFieldIn.setPtr(bFieldPtr, idx[0], idx[1]);
      bgKineticTempIn.setPtr(bgKineticTempPtr, idx[0], idx[1]);
      bgKineticNumDensityIn.setPtr(bgKineticNumDensityPtr, idx[0], idx[1]);
      perturbedKineticNumDensityIn.setPtr(perturbedKineticNumDensityPtr, idx[0], idx[1]);
      bgDistFIn.setPtr(bgDistFPtr, idx);
      perturbedFIn.setPtr(perturbedFPtr, idx);

      Eigen::VectorXd bFieldVec(nlocal2d);
      Eigen::VectorXd bgKineticTempVec(nlocal2d);
      Eigen::VectorXd bgKineticNumDensityVec(nlocal2d);
      Eigen::VectorXd perturbedKineticNumDensityVec(nlocal2d);
      Eigen::VectorXd bgDistFVec(nlocal4d);
      Eigen::VectorXd perturbedDistFVec(nlocal4d);

      for (int i = 0; i < nlocal2d; i++)
      {
        bFieldVec(i) = bFieldPtr[i];
        bgKineticTempVec(i) = bgKineticTempPtr[i];
        bgKineticNumDensityVec(i) = bgKineticNumDensityPtr[i];
        perturbedKineticNumDensityVec(i) = perturbedKineticNumDensityPtr[i];
      }

      for (int i = 0; i < nlocal4d; i++)
      {
        bgDistFVec(i) = bgDistFPtr[i];
        perturbedDistFVec(i) = perturbedFPtr[i];
      }

      // Interpolate data to quadrature points
      Eigen::VectorXd bFieldAtQuad = volQuad2d.interpMat*bFieldVec;
      Eigen::VectorXd bgKineticTempAtQuad = volQuad2d.interpMat*bgKineticTempVec;
      Eigen::VectorXd bgKineticNumDensityAtQuad = volQuad2d.interpMat*bgKineticNumDensityVec;
      Eigen::VectorXd perturbedKineticNumDensityAtQuad = volQuad2d.interpMat*perturbedKineticNumDensityVec;
      Eigen::VectorXd bgDistFAtQuad = volQuad4d.interpMat*bgDistFVec;
      Eigen::VectorXd perturbedDistFAtQuad = volQuad4d.interpMat*perturbedDistFVec;


      // perform quadrature
      for (int i = 0; i < nVolQuad4d; i++)
      {
        localInt += volQuad4d.weights[i]/(LX*LY)*( 2*Lucee::PI*bFieldAtQuad(i % nVolQuad2d)/kineticMass*
          bgKineticTempAtQuad(i % nVolQuad2d)*perturbedDistFAtQuad(i)*
          perturbedDistFAtQuad(i)/(2*bgDistFAtQuad(i)) +
          bgAdiabaticTemp/(2*bgKineticNumDensityAtQuad(i % nVolQuad2d))*
          perturbedKineticNumDensityAtQuad(i % nVolQuad2d)*
          perturbedKineticNumDensityAtQuad(i % nVolQuad2d) );
      }
    }

    double volInt = localInt;
    // get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    comm->allreduce(1, &localInt, &volInt, TX_SUM);
    
    std::vector<double> data(1);
    data[0] = volInt;
    fldInt.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  ETGFreeEnergy::declareTypes()
  {
    // Magnetic field B
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Background kinetic temperature
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Background kinetic number density
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Perturbed kinetic number density
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // Background kinetic distribution function
    this->appendInpVarType(typeid(Lucee::Field<4, double>));
    // Perturbed kinetic distribution function
    this->appendInpVarType(typeid(Lucee::Field<4, double>));

    // Output DynVector structure
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  ETGFreeEnergy::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
