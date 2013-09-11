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

    if (tbl.hasNumber("electronMass"))
      electronMass = tbl.getNumber("electronMass");
    else
      throw Lucee::Except("KineticTotalEnergyUpdater::readInput: Must specify electronMass");
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
    // zero-th moments of particle sources
    const Lucee::Field<1, double>& mom0SrcElcIn = this->getInp<Lucee::Field<1, double> >(4);
    const Lucee::Field<1, double>& mom0SrcIonIn = this->getInp<Lucee::Field<1, double> >(5);
    // second moments of particle sources
    const Lucee::Field<1, double>& mom2SrcElcIn = this->getInp<Lucee::Field<1, double> >(6);
    const Lucee::Field<1, double>& mom2SrcIonIn = this->getInp<Lucee::Field<1, double> >(7);
    // Electric potential
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(8);
    // Total heat fluxes at left and right edges
    const Lucee::DynVector<double>& heatFluxesIn = this->getInp<Lucee::DynVector<double> >(9);
    Lucee::DynVector<double>& totalEnergyOut = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();
    double dt = t - this->getCurrTime();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2ElcPtr = mom2ElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2IonPtr = mom2IonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom0SrcElcPtr = mom0SrcElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom0SrcIonPtr = mom0SrcIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2SrcElcPtr = mom2SrcElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2SrcIonPtr = mom2SrcIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    // Currently only using element[0], which is heat flux at right edge
    std::vector<double> heatFluxes = heatFluxesIn.getLastInsertedData();
    
    double totalEnergy = 0.0;
    double energyLost = ( heatFluxes[0] + heatFluxes[3] );
    double energyAdded = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      nElcIn.setPtr(nElcPtr, ix);
      nIonIn.setPtr(nIonPtr, ix);
      mom2ElcIn.setPtr(mom2ElcPtr, ix);
      mom2IonIn.setPtr(mom2IonPtr, ix);
      mom0SrcElcIn.setPtr(mom0SrcElcPtr, ix);
      mom0SrcIonIn.setPtr(mom0SrcIonPtr, ix);
      mom2SrcElcIn.setPtr(mom2SrcElcPtr, ix);
      mom2SrcIonIn.setPtr(mom2SrcIonPtr, ix);
      phiIn.setPtr(phiPtr, ix);

      Eigen::VectorXd totalEnergyAtNodes(nlocal);
      Eigen::VectorXd energyAddedAtNodes(nlocal);
      
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
      {
        // Compute the total energy at the element's nodes before interpolation
        double parallelEnergy = 0.5*electronMass*mom2ElcPtr[componentIndex] +
          0.5*ionMass*mom2IonPtr[componentIndex];
        double perpendicularEnergy = ELEMENTARY_CHARGE*tPerp*( nIonPtr[componentIndex] +
          nElcPtr[componentIndex] );
        double fieldEnergy = ELEMENTARY_CHARGE*phiPtr[componentIndex]*
          ( nIonPtr[componentIndex] - nElcPtr[componentIndex] );

        totalEnergyAtNodes(componentIndex) = parallelEnergy + perpendicularEnergy + fieldEnergy;

        // Compute the energy added to the system from the source at the element's nodes
        double parallelEnergyAdded = 0.5*electronMass*mom2SrcElcPtr[componentIndex] +
          0.5*ionMass*mom2SrcIonPtr[componentIndex];
        double perpendicularEnergyAdded = ELEMENTARY_CHARGE*tPerp*( mom0SrcElcPtr[componentIndex] +
            mom0SrcIonPtr[componentIndex] );
        double fieldEnergyAdded = ELEMENTARY_CHARGE*phiPtr[componentIndex]*
          ( mom0SrcIonPtr[componentIndex] - mom0SrcElcPtr[componentIndex] );

        energyAddedAtNodes(componentIndex) = ( parallelEnergyAdded + perpendicularEnergyAdded 
            + fieldEnergyAdded );
      }

      Eigen::VectorXd totalEnergyAtQuadPoints = interpMatrix*totalEnergyAtNodes;
      Eigen::VectorXd energyAddedAtQuadPoints = interpMatrix*energyAddedAtNodes;

      for (int componentIndex = 0; componentIndex < totalEnergyAtQuadPoints.rows(); componentIndex++)
      {
        totalEnergy += gaussWeights[componentIndex]*totalEnergyAtQuadPoints(componentIndex);
        energyAdded += gaussWeights[componentIndex]*energyAddedAtQuadPoints(componentIndex);
      }
    }

    std::vector<double> data(3);
    // Total energy at this time step
    data[0] = totalEnergy;
    // Energy lost at this time step
    data[1] = energyLost*dt;
    // Energy added at this time step
    data[2] = energyAdded*dt;

    // Put data into the DynVector
    totalEnergyOut.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  KineticTotalEnergyUpdater::declareTypes()
  {
    // takes inputs (n_e(x), n_i(x), <v^2>_e(x), <v^2>_i(x),
    // mom0srcElc, mom0srcIon, mom2srcElc, mom2srcIon, phi(x)
    // heat fluxes at edges
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
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
