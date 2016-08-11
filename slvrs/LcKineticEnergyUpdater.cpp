/**
 * @file	LcKineticEnergyUpdater.cpp
 *
 * @brief	Updater to compute energy using discrete hamiltonian and distribution function
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcKineticEnergyUpdater.h>
#include <LcMathPhysConstants.h>
// for output
#include <LcDynVector.h>

namespace Lucee
{
// set id for module system
  const char *KineticEnergyUpdater::id = "KineticEnergyUpdater";

  KineticEnergyUpdater::KineticEnergyUpdater()
    : UpdaterIfc()
  {
  }

  void 
  KineticEnergyUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<2> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<2> >("basis");
    else
      throw Lucee::Except("KineticEnergyUpdater::readInput: Must specify element to use using 'basis'");
  }

  void 
  KineticEnergyUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<2>& grid 
      = this->getGrid<Lucee::StructuredGridBase<2> >();
    // global region to update
    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<2> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[2];
    seq.fillWithIndex(idx);
    nodalBasis->setIndex(idx);
    
    int nlocal = nodalBasis->getNumNodes();
    int numQuadNodes = nodalBasis->getNumGaussNodes();

    // Get interpolation matrix, gaussian quadrature points, and weights
    Lucee::Matrix<double> interpMatrixLucee(numQuadNodes, nlocal);

    Lucee::Matrix<double> gaussOrdinatesLucee(numQuadNodes, 3);

    gaussWeights = std::vector<double>(numQuadNodes);

    // Allocate Eigen matrices
    interpMatrix = Eigen::MatrixXd(numQuadNodes, nlocal);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
  }

  Lucee::UpdaterStatus 
  KineticEnergyUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<2>& grid
      = this->getGrid<Lucee::StructuredGridBase<2> >();

    // Distribution function and hamiltonian
    const Lucee::Field<2, double>& distfIn = this->getInp<Lucee::Field<2, double> >(0);
    const Lucee::Field<2, double>& hamilIn = this->getInp<Lucee::Field<2, double> >(1);
    // Total kinetic energy
    Lucee::DynVector<double>& totalEnergyOut = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();
    double dt = t - this->getCurrTime();

    Lucee::Region<2, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> distfPtr = distfIn.createConstPtr();
    Lucee::ConstFieldPtr<double> hamilPtr = hamilIn.createConstPtr();

    double totalEnergy = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      for (int iv = globalRgn.getLower(1); iv < globalRgn.getUpper(1); iv++)
      {
        distfIn.setPtr(distfPtr, ix, iv);
        hamilIn.setPtr(hamilPtr, ix, iv);

        Eigen::VectorXd distfAtNodes(nlocal);
        Eigen::VectorXd hamilAtNodes(nlocal);

        // Fill out vectors at nodes
        for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        {
          distfAtNodes(componentIndex) = distfPtr[componentIndex];
          hamilAtNodes(componentIndex) = hamilPtr[componentIndex];
        }

        // Interpolate both vectors to gaussian quadrature points
        Eigen::VectorXd distfAtQuadPoints = interpMatrix*distfAtNodes;
        Eigen::VectorXd hamilAtQuadPoints = interpMatrix*hamilAtNodes;

        // Integrate hamil*distf over this cell and add to total
        for (int quadIndex = 0; quadIndex < distfAtQuadPoints.rows(); quadIndex++)
          totalEnergy += gaussWeights[quadIndex]*distfAtQuadPoints(quadIndex)*hamilAtQuadPoints(quadIndex);
      }
    }

    std::vector<double> data(1);
    data[0] = totalEnergy;

    // Put data into the DynVector
    totalEnergyOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  KineticEnergyUpdater::declareTypes()
  {
    // Takes distf and hamil as inputs
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    this->appendInpVarType(typeid(Lucee::Field<2, double>));
    // returns one output: dynvector containing kinetic energy at time t
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  KineticEnergyUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
