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

    if (tbl.hasNumber("B0"))
      B0 = tbl.getNumber("B0");
    else
      throw Lucee::Except("KineticHeatFluxAtEdge3DUpdater::readInput: Must specify B0 using 'B0'");

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
    // Returns heat flux vs time as a dynvector
    Lucee::DynVector<double>& qVsTime = this->getOut<Lucee::DynVector<double> >(0);

    int nlocal = nodalBasis->getNumNodes();

    Lucee::Region<1, int> globalRgn = grid.getGlobalRegion();

    Lucee::ConstFieldPtr<double> phiPtr = phiIn.createConstPtr();

    std::vector<double> momentsAtEdgesElc = momentsAtEdgesElcIn.getLastInsertedData();
    std::vector<double> momentsAtEdgesIon = momentsAtEdgesIonIn.getLastInsertedData();
    
    // Find value of the following input fields at the right-most edge of the domain
    phiIn.setPtr(phiPtr, globalRgn.getUpper(0)-1);

    // Guide:
    // momentsAtEdges[0 to 2] = <vPara>,<vPara^3>,<mu*vPara> moments of f on right edge only
    
    double ionHeatFluxRight = 0.5*ionMass*momentsAtEdgesIon[1] + 
      momentsAtEdgesIon[0]*ELEMENTARY_CHARGE*phiPtr[nlocal-1] +
      momentsAtEdgesIon[2]*B0;

    double elcHeatFluxRight = 0.5*elcMass*momentsAtEdgesElc[1] - 
      momentsAtEdgesElc[0]*ELEMENTARY_CHARGE*phiPtr[nlocal-1] +
      momentsAtEdgesElc[2]*B0;

    /*if (elcHeatFluxRight < 0.0)
    {
      std::cout << "term1 = " << 0.5*elcMass*momentsAtEdgesElc[1] << std::endl;
      std::cout << "term2 = " << -momentsAtEdgesElc[0]*ELEMENTARY_CHARGE*phiPtr[nlocal-1] << std::endl;
      std::cout << "term3 = " << momentsAtEdgesElc[2]*B0 << std::endl;
    }*/

    std::vector<double> data(3);
    data[0] = ionHeatFluxRight + elcHeatFluxRight;
    data[1] = ionHeatFluxRight;
    data[2] = elcHeatFluxRight;

    //std::cout << "momentsAtEdgesElc[0] = " << momentsAtEdgesElc[0] << std::endl;
    //std::cout << "momentsAtEdgesIon[0] = " << momentsAtEdgesIon[0] << std::endl;
    
    qVsTime.appendData(t, data);

    if (computeSheathCoefficient == true)
    {
      // tPerp*n fields of both species
      const Lucee::Field<1, double>& tElcIn = this->getInp<Lucee::Field<1, double> >(3);
      const Lucee::Field<1, double>& tIonIn = this->getInp<Lucee::Field<1, double> >(4);
      Lucee::DynVector<double>& sheathCoefficientVsTime = this->getOut<Lucee::DynVector<double> >(1);
       
      Lucee::ConstFieldPtr<double> tElcPtr = tElcIn.createConstPtr();
      Lucee::ConstFieldPtr<double> tIonPtr = tIonIn.createConstPtr();
    
      tElcIn.setPtr(tElcPtr, globalRgn.getUpper(0)-1);
      tIonIn.setPtr(tIonPtr, globalRgn.getUpper(0)-1);
      
      std::vector<double> sheathData(5);
      
      double ionHeatFluxPresheath = 0.5*ionMass*momentsAtEdgesIon[1] + 
        momentsAtEdgesIon[2]*B0;

      double elcHeatFluxPresheath = 0.5*elcMass*momentsAtEdgesElc[1] + 
        momentsAtEdgesElc[2]*B0;
      double kTe0 = tElcPtr[nlocal-1]*ELEMENTARY_CHARGE;
      double kTi0 = tIonPtr[nlocal-1]*ELEMENTARY_CHARGE;
      // Total transmission coefficient
      sheathData[0] = (ionHeatFluxPresheath+elcHeatFluxPresheath)/( kTe0*momentsAtEdgesElc[0] );
      // Ion transmission coefficient
      sheathData[1] = ionHeatFluxPresheath/( kTi0*momentsAtEdgesIon[0] );
      // Elc transmission coefficient
      sheathData[2] = elcHeatFluxPresheath/( kTe0*momentsAtEdgesElc[0] );
      sheathData[3] = tElcPtr[nlocal-1];
      sheathData[4] = tIonPtr[nlocal-1];
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
    // inputs: elc and ion temperatures
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    this->appendInpVarType(typeid(Lucee::Field<1, double>));
    // returns one output: dynvector of heat flux at edge vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
    // optional second output: dynvector of sheath power transmission
    // coefficient vs time
    this->appendOutVarType(typeid(Lucee::DynVector<double>));
  }
}
