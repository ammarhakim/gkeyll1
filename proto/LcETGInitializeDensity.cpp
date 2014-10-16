/**
 * @file	LcETGInitializeDensity.cpp
 *
 * @brief	Updater to compute 2d moments of a 4d distribution function with an additional weighting function.
 * Currently only works for 0th nIn
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcETGInitializeDensity.h>

// std includes
#include <vector>

namespace Lucee
{
  const char *ETGInitializeDensity::id = "ETGInitializeDensity";

  ETGInitializeDensity::ETGInitializeDensity()
    : Lucee::UpdaterIfc()
  {
  }
  
  void
  ETGInitializeDensity::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasNumber("polyOrder"))
      polyOrder = tbl.getNumber("polyOrder");
    else
      throw Lucee::Except("ETGInitializeDensity::readInput: Must specify polyOrder");

    // get hold of 4D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<4> >("basis4d"))
      nodalBasis4d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<4> >("basis4d");
    else
      throw Lucee::Except(
        "ETGInitializeDensity::readInput: Must specify 4D element to use using 'basis4d'");

    // get hold of 2D element to use
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis2d"))
      nodalBasis2d = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis2d");
    else
      throw Lucee::Except(
        "ETGInitializeDensity::readInput: Must specify 2D element to use using 'basis1d'");

    // get nIn to compute
    if (tbl.hasNumber("constantDensity"))
    constantDensity = tbl.getNumber("constantDensity");
    else
      throw Lucee::Except(
        "ETGInitializeDensity::readInput: Must specify density using 'constantDensity'");
  }

  void
  ETGInitializeDensity::initialize()
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

    if (polyOrder == 1)
    {
      mappingMatrix = Eigen::MatrixXd(16,4);
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
    }
    else if (polyOrder == 2)
    {
      mappingMatrix = Eigen::MatrixXd(48,8);
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
                       0,0,0,0,0,0,0,1,
                       1,0,0,0,0,0,0,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,0,1,0,0,0,
                       0,0,0,0,0,0,1,0,
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
  }

  Lucee::UpdaterStatus
  ETGInitializeDensity::update(double t)
  {
    
    // get hold of grid
    const Lucee::StructuredGridBase<4>& grid 
      = this->getGrid<Lucee::StructuredGridBase<4> >();

    // get input field (2d)
    const Lucee::Field<2, double>& nIn = this->getInp<Lucee::Field<2, double> >(0);
    // get output field (2d)
    Lucee::Field<4, double>& distF = this->getOut<Lucee::Field<4, double> >(0);

    // local region to update (This is the 4D region. The 2D region is
    // assumed to have the same cell layout as the X-direction of the 4D region)
    Lucee::Region<4, int> localExtRgn = distF.getExtRegion();

    // iterators into fields
    Lucee::ConstFieldPtr<double> nPtr = nIn.createConstPtr();
    Lucee::FieldPtr<double> distFPtr = distF.createPtr();

    int idx[4];
    Lucee::RowMajorSequencer<4> seq(localExtRgn);
    unsigned nlocal2d = nodalBasis2d->getNumNodes();
    unsigned nlocal4d = nodalBasis4d->getNumNodes();
    int nVolQuad2d = nodalBasis2d->getNumGaussNodes();
    // Loop over each cell in 4D space
    while(seq.step())
    {
      seq.fillWithIndex(idx);

      grid.setIndex(idx);

      nIn.setPtr(nPtr, idx[0], idx[1]);
      distF.setPtr(distFPtr, idx);

      Eigen::VectorXd scaleFactors2d(nlocal2d);
      for (int i = 0; i < nlocal2d; i++)
        scaleFactors2d(i) = constantDensity/nPtr[i];

      Eigen::VectorXd scaleFactors4d = mappingMatrix*scaleFactors2d;

      for (int i = 0; i < nlocal4d; i++)
        distFPtr[i] = scaleFactors4d(i)*distFPtr[i];
    }

    return Lucee::UpdaterStatus();
  }

  void
  ETGInitializeDensity::declareTypes()
  {
    // computed density (2d)
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // output scaled distribution function (4d)
    this->appendOutVarType(typeid(Lucee::Field<4, double>));
  }

  void
  ETGInitializeDensity::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
