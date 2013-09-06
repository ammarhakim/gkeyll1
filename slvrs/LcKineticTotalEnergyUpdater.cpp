/**
 * @file	LcKineticTotalEnergyUpdater.cpp
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcField.h>
#include <LcLinAlgebra.h>
#include <LcMathLib.h>
#include <LcKineticTotalEnergyUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// for cutoff velocities
#include <LcDynVector.h>

// math include
#include <cmath>

namespace Lucee
{
// set id for module system
  const char *KineticTotalEnergyUpdater::id = "KineticTotalEnergyUpdater";

  KineticTotalEnergyUpdater::KineticTotalEnergyUpdater()
    : UpdaterIfc()
  {
  }

  void 
  KineticTotalEnergyUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("tPerp"))
      tPerp = tbl.getNumber("tPerp");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify tPerp.");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify ionMass");
  }

  void 
  KineticTotalEnergyUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();

    // get hold of grid
    const Lucee::StructuredGridBase<1>& grid 
      = this->getGrid<Lucee::StructuredGridBase<1> >();
    // global region to update
    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::RowMajorSequencer<1> seq(globalRgn);
    seq.step(); // just to get to first index
    int idx[1];
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
    gaussOrdinates = Eigen::MatrixXd(numQuadNodes, 3);

    // Get the interpolation matrix for the volume quadrature points
    nodalBasis->getGaussQuadData(interpMatrixLucee, gaussOrdinatesLucee, gaussWeights);

    copyLuceeToEigen(interpMatrixLucee, interpMatrix);
    copyLuceeToEigen(gaussOrdinatesLucee, gaussOrdinates);
  }

  Lucee::UpdaterStatus 
  KineticTotalEnergyUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // Electron and ion densities
    const Lucee::Field<1, double>& nElcIn = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(1);
    // Electron and ion moment fluxes
    const Lucee::Field<1, double>& mom2ElcIn = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& mom2IonIn = this->getInp<Lucee::Field<1, double> >(3);
    // Electric potential
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(4);
    // Total heat fluxes at left and right edges
    const Lucee::DynVector<double>& heatFluxesIn = this->getInp<Lucee::DynVector<double> >(5);
    Lucee::DynVector<double>& totalEnergyOut = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2ElcPtr = mom2ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2IonPtr = mom2IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    // Currently only using element[0], which is heat flux at right edge
    std::vector<double> heatFluxes = heatFluxesIn.getLastInsertedData();
    
    double totalEnergy = 0.0;
    double energyLost = 0.0;
    double energyAdded = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nElcIn.setPtr(nElcPtr, ix);
      nIonIn.setPtr(nIonPtr, ix);
      mom2ElcIn.setPtr(mom2ElcPtr, ix);
      mom2IonIn.setPtr(mom2IonPtr, ix);
      phiIn.setPtr(phiPtr, ix);

      Eigen::VectorXd totalEnergyAtNodes(nlocal);
      
      // Compute the total energy at the element's nodes before interpolation
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        double parallelEnergy = 0.5*ELECTRON_MASS*mom2ElcPtr[componentIndex] +
          0.5*ionMass*mom2IonPtr[componentIndex];
        double perpendicularEnergy = 4*ELEMENTARY_CHARGE*tPerp*(nIonPtr[componentIndex] +
          nElcPtr[componentIndex]);
        double fieldEnergy = ELEMENTARY_CHARGE*(nIonPtr[componentIndex] - nElcPtr[componentIndex])*
          phiPtr[componentIndex];

        totalEnergyAtNodes(componentIndex) = parallelEnergy + perpendicularEnergy + fieldEnergy;
      }

      Eigen::VectorXd totalEnergyAtQuadPoints = interpMatrix*totalEnergyAtNodes;

      for (int componentIndex = 0; componentIndex < totalEnergyAtQuadPoints.rows(); componentIndex++)
        totalEnergy += gaussWeights[componentIndex]*totalEnergyAtQuadPoints(componentIndex);
    }

    std::vector<double> data(3);
    data[0] = totalEnergy;
    data[1] = energyLost;
    data[2] = energyAdded;

    // Put data into the DynVector
    totalEnergyOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  KineticTotalEnergyUpdater::declareTypes()
  {
    // takes inputs (n_e(x), n_i(x), <v^2>_e(x), <v^2>_i(x), heat fluxes at edges
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // returns one output: dynvector containing energy at time t
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }

  void
  KineticTotalEnergyUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
