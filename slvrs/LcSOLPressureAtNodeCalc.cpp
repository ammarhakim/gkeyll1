/**
 * @file	LcSOLPressureAtNodeCalc.cpp
 *
 * @brief	Simple updater to compute pressure node-by-node for 5D SOL simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLPressureAtNodeCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLPressureAtNodeCalc::id = "SOLPressureAtNodeCalc";

  SOLPressureAtNodeCalc::SOLPressureAtNodeCalc()
  {
  }

  void
  SOLPressureAtNodeCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOLPressureAtNodeCalc::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("SOLPressureAtNodeCalc::readInput: Must specify electron mass using 'speciesMass'");
  }

  void
  SOLPressureAtNodeCalc::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLPressureAtNodeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    const Lucee::Field<3, double>& numDensIn = this->getInp<Lucee::Field<3, double> >(0);
    const Lucee::Field<3, double>& energyIn = this->getInp<Lucee::Field<3, double> >(1);
    const Lucee::Field<3, double>& meanVelocityIn = this->getInp<Lucee::Field<3, double> >(2);
    Lucee::Field<3, double>& pressureOut = this->getOut<Lucee::Field<3, double> >(0);

    Lucee::ConstFieldPtr<double> numDensPtr = numDensIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyPtr = energyIn.createConstPtr();
    Lucee::ConstFieldPtr<double> meanVelocityPtr = meanVelocityIn.createConstPtr();
    
    Lucee::FieldPtr<double> pressurePtr = pressureOut.createPtr();

    unsigned nlocal = nodalBasis->getNumNodes();

    double cellCentroid[3];
    int idx[3];
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<3> seq(localRgn);
    bool negativePressureStatus = false;

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      
      numDensIn.setPtr(numDensPtr, idx);
      energyIn.setPtr(energyPtr, idx);
      meanVelocityIn.setPtr(meanVelocityPtr, idx);
      
      pressureOut.setPtr(pressurePtr, idx);

      for (int i = 0; i < nlocal; i++)
        pressurePtr[i] = 2.0/3.0*(energyPtr[i] - 
          0.5*speciesMass*meanVelocityPtr[i]*meanVelocityPtr[i]/numDensPtr[i]);
    }

    return Lucee::UpdaterStatus();
  }

  void
  SOLPressureAtNodeCalc::declareTypes()
  {
    // Input: nodal number density
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: nodal energy
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: nodal <v> = nu
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: nodal pressure
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
