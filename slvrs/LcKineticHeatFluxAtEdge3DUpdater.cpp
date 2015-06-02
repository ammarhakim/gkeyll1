/**
 * @file	LcKineticKineticHeatFluxAtEdge3DUpdater.cpp
 *
 * @brief	Updater to compute heat flux at right-most edge in domain
 * Used for kinetic SOL problem with 1D2V
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcKineticHeatFluxAtEdge3DUpdater.h>
#include <LcMathPhysConstants.h>
// For function handle:
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for module system
  const char *KineticHeatFluxAtEdge3DUpdater::id = "KineticHeatFluxAtEdge3DUpdater";

  KineticHeatFluxAtEdge3DUpdater::KineticHeatFluxAtEdge3DUpdater()
    : UpdaterIfc()
  {
  }

  void 
  KineticHeatFluxAtEdge3DUpdater::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<1> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<1> >("basis");
    else
      throw Lucee::Except("KineticHeatFluxAtEdge3DUpdater::readInput: Must specify element to use using 'basis1d'");

    if (tbl.hasNumber("ionMass"))
      ionMass = tbl.getNumber("ionMass");
    else
      throw Lucee::Except("KineticHeatFluxAtEdge3DUpdater::readInput: Must specify ionMass");

    if (tbl.hasNumber("electronMass"))
      elcMass = tbl.getNumber("electronMass");
    else
      throw Lucee::Except("KineticHeatFluxAtEdge3DUpdater::readInput: Must specify electronMass");

    // Check to see if we need to compute the sheat power transmission coefficient
    computeSheathCoefficient = false;
    if (tbl.hasBool("computeSheathCoefficient"))
      computeSheathCoefficient = tbl.getBool("computeSheathCoefficient");
  }

  void 
  KineticHeatFluxAtEdge3DUpdater::initialize()
  {
    Lucee::UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus 
  KineticHeatFluxAtEdge3DUpdater::update(double t)
  {
    const Lucee::StructuredGridBase<1>& grid
      = this->getGrid<Lucee::StructuredGridBase<1> >();

    // phi
    const Lucee::Field<1, double>& phiIn = this->getInp<Lucee::Field<1, double> >(0);
    // parallel velocity moments at left and right edges
    const Lucee::DynVector<double>& momentsAtEdgesElcIn = this->getInp<Lucee::DynVector<double> >(1);
    const Lucee::DynVector<double>& momentsAtEdgesIonIn = this->getInp<Lucee::DynVector<double> >(2);
    // number densities of both species
    const Lucee::Field<1, double>& nIonIn = this->getInp<Lucee::Field<1, double> >(3);
    const Lucee::Field<1, double>& nElcIn = this->getInp<Lucee::Field<1, double> >(4);
    // tPerp*n fields of both species
    const Lucee::Field<1, double>& tPerpElcIn = this->getInp<Lucee::Field<1, double> >(5);
    const Lucee::Field<1, double>& tPerpIonIn = this->getInp<Lucee::Field<1, double> >(6);
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nElcPtr = nElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> nIonPtr = nIonIn.createConstPtr();
    Lucee::ConstFieldPtr<double> tPerpElcPtr = tPerpElcIn.createConstPtr();
    Lucee::ConstFieldPtr<double> tPerpIonPtr = tPerpIonIn.createConstPtr();

    std::vector<double> momentsAtEdgesElc = momentsAtEdgesElcIn.getLastInsertedData();
    std::vector<double> momentsAtEdgesIon = momentsAtEdgesIonIn.getLastInsertedData();
    
    // Find value of the following input fields at the right-most edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);
    nElcIn.setPtr(nElcPtr, globalRgn.getUpper(0)-1);
    nIonIn.setPtr(nIonPtr, globalRgn.getUpper(0)-1);
    tPerpElcIn.setPtr(tPerpElcPtr, globalRgn.getUpper(0)-1);
    tPerpIonIn.setPtr(tPerpIonPtr, globalRgn.getUpper(0)-1);

    // Perpendicular electron and ion temperatures (in eV)
    double tPerpElc = tPerpElcPtr[nlocal-1]/nElcPtr[nlocal-1];
    double tPerpIon = tPerpIonPtr[nlocal-1]/nIonPtr[nlocal-1];

    // Guide:
    // momentsAtEdges[0 to 3] mom 0 to 3 on right edge
    // momentsAtEdges[4 to 7] mom 0 to 3 on left edge
    
    double ionParallelHeatFluxRight = 0.5*ionMass*momentsAtEdgesIon[7] + 
      momentsAtEdgesIon[5]*ELEMENTARY_CHARGE*phiPtr[nlocal-1];
    double ionHeatFluxRight = ionParallelHeatFluxRight +
      momentsAtEdgesIon[5]*ELEMENTARY_CHARGE*tPerpIon;

    double elcParallelHeatFluxRight = 0.5*elcMass*momentsAtEdgesElc[7] - 
      momentsAtEdgesElc[5]*ELEMENTARY_CHARGE*phiPtr[nlocal-1];
    double elcHeatFluxRight = elcParallelHeatFluxRight + 
      momentsAtEdgesElc[5]*ELEMENTARY_CHARGE*tPerpElc;

    // Find value of the following input fields at the left-most edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getLower(0));

    double ionParallelHeatFluxLeft = 0.5*ionMass*momentsAtEdgesIon[3] + 
      momentsAtEdgesIon[1]*ELEMENTARY_CHARGE*phiPtr[0];
    double ionHeatFluxLeft = ionParallelHeatFluxLeft +
      momentsAtEdgesIon[1]*ELEMENTARY_CHARGE*tPerpIon;

    double elcParallelHeatFluxLeft = 0.5*elcMass*momentsAtEdgesElc[3] - 
      momentsAtEdgesElc[1]*ELEMENTARY_CHARGE*phiPtr[0];
    double elcHeatFluxLeft = elcParallelHeatFluxLeft + 
      momentsAtEdgesElc[1]*ELEMENTARY_CHARGE*tPerpElc;

    std::vector<double> data(6);
    data[0] = ionHeatFluxRight + elcHeatFluxRight;
    data[1] = ionHeatFluxRight;
    data[2] = elcHeatFluxRight;
    data[3] = -(ionHeatFluxLeft + elcHeatFluxLeft);
    data[4] = -ionHeatFluxLeft;
    data[5] = -elcHeatFluxLeft;
    
    qVsTime.appendData(t, data);

    if (computeSheathCoefficient == true)
    {
      Lucee::DynVector<double>& sheathCoefficientVsTime = this->getOut<Lucee::DynVector<double> >(1);
      
      std::vector<double> sheathData(7);

      double kTe0 = elcMass*(momentsAtEdgesElc[6] -
          momentsAtEdgesElc[5]*momentsAtEdgesElc[5]/momentsAtEdgesElc[4])/
          momentsAtEdgesElc[4];
      
      double kTi0 = ionMass*(momentsAtEdgesIon[6] - 
          momentsAtEdgesIon[5]*momentsAtEdgesIon[5]/momentsAtEdgesIon[4])/
          momentsAtEdgesIon[4];

      // Total transmission coefficient
      sheathData[0] = data[0]/( kTe0*momentsAtEdgesElc[5] );
      sheathData[1] = data[1]/( kTi0*momentsAtEdgesIon[5] );
      sheathData[2] = data[2]/( kTe0*momentsAtEdgesElc[5] );
      sheathData[3] = kTi0;
      sheathData[4] = kTe0;
      // Diagnostics
      sheathData[5] = momentsAtEdgesElc[5];
      sheathData[6] = momentsAtEdgesIon[5];

      sheathCoefficientVsTime.appendData(t, sheathData);
    }

    return Lucee::UpdaterStatus();
  }

  void
  KineticHeatFluxAtEdge3DUpdater::declareTypes()
  {
    // input: phi(x)
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // inputs: elc and ion dynvectors containing moments at edges of domain
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    this->appendInpVarType(typeid(Lucee::DynVector<double>));
    // inputs: elc and ion number densitites
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // inputs: elc and ion perpendicular temperature * number densities
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
    // optional second output: dynvector of sheath power transmission
    // coefficient vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
}
