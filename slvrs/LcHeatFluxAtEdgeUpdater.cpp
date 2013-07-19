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

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("HeatFluxAtEdgeUpdater::readInput: Must specify speciesMass");
  }

  void 
  HeatFluxAtEdgeUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  HeatFluxAtEdgeUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // First three velocity moments + driftU
    const Lucee::Field<1, double>& mom1In = this->getInp<Lucee::Field<1, double> >(0);
    const Lucee::Field<1, double>& mom2In = this->getInp<Lucee::Field<1, double> >(1);
    const Lucee::Field<1, double>& mom3In = this->getInp<Lucee::Field<1, double> >(2);
    const Lucee::Field<1, double>& driftUIn = this->getInp<Lucee::Field<1, double> >(3);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> mom1Ptr = mom1In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2Ptr = mom2In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom3Ptr = mom3In.createConstPtr();
    Lucee::ConstFieldPtr<double> driftUPtr = driftUIn.createConstPtr();

    // Find value of the following input fields at the very last edge of the domain
    mom1In.setPtr(mom1Ptr, globalRgn.getUpper(0)-1);
    mom2In.setPtr(mom2Ptr, globalRgn.getUpper(0)-1);
    mom3In.setPtr(mom3Ptr, globalRgn.getUpper(0)-1);
    driftUIn.setPtr(driftUPtr, globalRgn.getUpper(0)-1);
    
    double heatFlux = 0.5*speciesMass*(mom3Ptr[nlocal-1] 
      - 3*driftUPtr[nlocal-1]*mom2Ptr[nlocal-1]
      + 3*driftUPtr[nlocal-1]*driftUPtr[nlocal-1]*mom1Ptr[nlocal-1]
      - driftUPtr[nlocal-1]*driftUPtr[nlocal-1]*driftUPtr[nlocal-1]);
    
    std::vector<double> data(1);
    data[0] = heatFlux;
    qVsTime.appendData(t, data);

    return Lucee::UpdaterStatus();
  }

  void
  HeatFluxAtEdgeUpdater::declareTypes()
  {
    // takes four inputs (<v>,<v^2>,<v^3>) moments + driftU
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
}
