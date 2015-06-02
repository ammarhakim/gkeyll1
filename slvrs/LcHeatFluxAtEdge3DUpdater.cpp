/**
 * @file	LcHeatFluxAtEdge3DUpdater.cpp
 *
 * @brief	Updater to compute heat flux at right-most edge in domain
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
#include <LcHeatFluxAtEdge3DUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>
// For function handle:
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *HeatFluxAtEdge3DUpdater::id = "HeatFluxAtEdge3DUpdater";

  HeatFluxAtEdge3DUpdater::HeatFluxAtEdge3DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  HeatFluxAtEdge3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("HeatFluxAtEdge3DUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("HeatFluxAtEdge3DUpdater::readInput: Must specify ionMass");
  }

  void 
  HeatFluxAtEdge3DUpdater::initialize()
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

    // Get mass matrix and then copy to Eigen format
    Lucee::Matrix<double> massMatrixLucee(nlocal, nlocal);

    nodalBasis->getMassMatrix(massMatrixLucee);
    
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
  HeatFluxAtEdge3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // phi
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(0);
    // parallel velocity moments at left and right edges
    const Lucee::DynVector<double>& momentsAtEdgesIn = this->getInp<Lucee::DynVector<double> >(1);
    // vThermSq integrated in domain
    const Lucee::DynVector<double>& integratedVtSqIn = this->getInp<Lucee::DynVector<double> >(2);
    // nIon field
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(3);
    // tPerp field
    const Lucee::Field<1, double>& tPerpIn = this->getInp<Lucee::Field<1, double> >(4);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> tPerpPtr = tPerpIn.createConstPtr();

    std::vector<double> momentsAtEdges = momentsAtEdgesIn.getLastInsertedData();
    std::vector<double> integratedVtSq = integratedVtSqIn.getLastInsertedData();
    
    // Divide integrated vtSq (to get T_parallel) by length of domain
    // Consider using grid.getNumCells(0)
    double meanVtSq = integratedVtSq[0]/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

    // Find value of the following input fields at the very last edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);
    nIonIn.setPtr(nIonPtr, globalRgn.getUpper(0)-1);
    tPerpIn.setPtr(tPerpPtr, globalRgn.getUpper(0)-1);

    double tPerpIon = tPerpPtr[nlocal-1]/nIonPtr[nlocal-1];
    double tPerpElc = tPerpPtr[nlocal-1]/nIonPtr[nlocal-1];

    // Figure out value of perpendicular electron and ion temperatures
    double ionHeatFlux = 0.5*ionMass*momentsAtEdges[7] + momentsAtEdges[5]*ELEMENTARY_CHARGE*(tPerpIon + phiPtr[nlocal-1]);
    // Te = average of Ti(x)
    double electronHeatFlux = (ionMass*meanVtSq + ELEMENTARY_CHARGE*tPerpElc)*momentsAtEdges[5];
    
    std::vector<double> data(3);
    data[0] = ionHeatFlux + electronHeatFlux;
    data[1] = ionHeatFlux;
    data[2] = electronHeatFlux;
    qVsTime.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  HeatFluxAtEdge3DUpdater::declareTypes()
  {
    // input: phi(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // input: moments evaluated at edges input
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // input: vThermSq integrated in entire domain (need to divide by length in code)
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // input: ion density
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // input: tPerp (in eV)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
  
  void
  HeatFluxAtEdge3DUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
