/**
 * @file	LcHeatFluxAtEdgeUpdater.cpp
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
#include <LcHeatFluxAtEdgeUpdater.h>
#include <LcStructuredGridBase.h>
#include <LcMathPhysConstants.h>

namespace Lucee
{
// set id for module system
  const char *HeatFluxAtEdgeUpdater::id = "HeatFluxAtEdgeUpdater";

  HeatFluxAtEdgeUpdater::HeatFluxAtEdgeUpdater()
    : UpdaterIfc()
  {
  }

  void 
  HeatFluxAtEdgeUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("tPerp"))
      tPerp = tbl.getNumber("tPerp");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify tPerp");
  }

  void 
  HeatFluxAtEdgeUpdater::initialize()
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
  HeatFluxAtEdgeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // First three velocity moments + driftU
    const Lucee::Field<1, double>& mom1In = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& mom3In = this->getInp<Lucee::Field<1, double> >(1);
    const Lucee::Field<1, double>& vtSqIn = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(3);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> mom1Ptr = mom1In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom3Ptr = mom3In.createConstPtr();
    Lucee::ConstFieldPtr<double> vtSqPtr = vtSqIn.createConstPtr();
    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();
    
    // Find mean vtSq in entire domain
    double meanVtSq = 0.0;

    for (int ix = globalRgn.getLower(0); ix < globalRgn.getUpper(0); ix++)
    {
      // Set inputs
      vtSqIn.setPtr(vtSqPtr, ix);
      Eigen::VectorXd vtSqVec(nlocal);
      
      // Figure out vtSq(x) at quadrature points in the cell
      for (int componentIndex = 0; componentIndex < nlocal; componentIndex++)
        vtSqVec(componentIndex) = vtSqPtr[componentIndex];
      Eigen::VectorXd vtSqAtQuadPoints = interpMatrix*vtSqVec;

      for (int componentIndex = 0; componentIndex < vtSqAtQuadPoints.rows(); componentIndex++)
        meanVtSq += gaussWeights[componentIndex]*vtSqAtQuadPoints(componentIndex);
    }

    // Divide by length of domain
    // Consider using grid.getNumCells(0)
    meanVtSq = meanVtSq/(grid.getDx(0)*(globalRgn.getUpper(0)-globalRgn.getLower(0)));

    // Find value of the following input fields at the very last edge of the domain
    mom1In.setPtr(mom1Ptr, globalRgn.getUpper(0)-1);
    mom3In.setPtr(mom3Ptr, globalRgn.getUpper(0)-1);
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);
    
    double ionHeatFlux = 0.5*ionMass*mom3Ptr[nlocal-1] + mom1Ptr[nlocal-1]*ELEMENTARY_CHARGE*(tPerp + phiPtr[nlocal-1]);
    double electronHeatFlux = (ionMass*meanVtSq + ELEMENTARY_CHARGE*tPerp)*mom1Ptr[nlocal-1];
    
    std::vector<double> data(3);
    data[0] = ionHeatFlux + electronHeatFlux;
    data[1] = ionHeatFlux;
    data[2] = electronHeatFlux;
    qVsTime.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  HeatFluxAtEdgeUpdater::declareTypes()
  {
    // takes four inputs (<v>,<v^3>) moments + vtSq + phi(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
  
  void
  HeatFluxAtEdgeUpdater::copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix,
    Eigen::MatrixXd& destinationMatrix)
  {
    for (int rowIndex = 0; rowIndex < destinationMatrix.rows(); rowIndex++)
      for (int colIndex = 0; colIndex < destinationMatrix.cols(); colIndex++)
        destinationMatrix(rowIndex, colIndex) = sourceMatrix(rowIndex, colIndex);
  }
}
