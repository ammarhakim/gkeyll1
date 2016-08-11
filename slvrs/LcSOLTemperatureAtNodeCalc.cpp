/**
 * @file	LcSOLTemperatureAtNodeCalc.cpp
 *
 * @brief	Simple updater to compute temperature node-by-node for 5D SOL simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLTemperatureAtNodeCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLTemperatureAtNodeCalc::id = "SOLTemperatureAtNodeCalc";

  SOLTemperatureAtNodeCalc::SOLTemperatureAtNodeCalc()
  {
  }

  void
  SOLTemperatureAtNodeCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOLTemperatureAtNodeCalc::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("SOLTemperatureAtNodeCalc::readInput: Must specify electron mass using 'speciesMass'");
  }

  void
  SOLTemperatureAtNodeCalc::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLTemperatureAtNodeCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Input number density of species weighted by 2*pi/m
    const Lucee::Field<3, double>& numDensIn = this->getInp<Lucee::Field<3, double> >(0);
    // Input: <H*f> weighted by 2*pi/m
    const Lucee::Field<3, double>& energyIn = this->getInp<Lucee::Field<3, double> >(1);
    // Input: <dH/dv*f> weighted by 2*pi/m
    const Lucee::Field<3, double>& meanVelocityIn = this->getInp<Lucee::Field<3, double> >(2);
    // Output dynvector of electron temperatures at wall
    Lucee::Field<3, double>& temperatureOut = this->getOut<Lucee::Field<3, double> >(0);

    Lucee::ConstFieldPtr<double> numDensPtr = numDensIn.createConstPtr();
    Lucee::ConstFieldPtr<double> energyPtr = energyIn.createConstPtr();
    Lucee::ConstFieldPtr<double> meanVelocityPtr = meanVelocityIn.createConstPtr();
    
    Lucee::FieldPtr<double> temperaturePtr = temperatureOut.createPtr();

    unsigned nlocal = nodalBasis->getNumNodes();

    double cellCentroid[3];
    int idx[3];
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<3> seq(localRgn);
    bool negativeTemperatureStatus = false;

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      
      numDensIn.setPtr(numDensPtr, idx);
      energyIn.setPtr(energyPtr, idx);
      meanVelocityIn.setPtr(meanVelocityPtr, idx);
      
      temperatureOut.setPtr(temperaturePtr, idx);

      for (int i = 0; i < nlocal; i++)
      {
        temperaturePtr[i] = 2.0/3.0*(energyPtr[i]/numDensPtr[i] - 
          0.5*speciesMass*meanVelocityPtr[i]*meanVelocityPtr[i]/(numDensPtr[i]*numDensPtr[i]));
        if (temperaturePtr[i] < 0.0)
        {
          negativeTemperatureStatus = true;
          // Print out information
          std::cout << "temperaturePtr[" << i << "] = " << temperaturePtr[i] << std::endl;
          std::cout << "idx = (" << idx[0] << "," << idx[1] << "," << idx[2] << ")" << std::endl;
          std::cout << "n = " << numDensPtr[i] << std::endl;
          std::cout << "E = " << energyPtr[i] << std::endl;
          std::cout << "u = " << meanVelocityPtr[i] << std::endl;
          std::cout << "m = " << speciesMass << std::endl;
        }
        else if (std::isnan(temperaturePtr[i]))
        {
          //std::cout << "Temperature is zero at node " << i << " in cell "
          //  << idx[0] << "," << idx[1] << "," << idx[2] << std::endl;
          temperaturePtr[i] = 0.0;
        }
      }
    }

    if (negativeTemperatureStatus == true)
      return Lucee::UpdaterStatus(false, 0.0);
    else
      return Lucee::UpdaterStatus();
  }

  void
  SOLTemperatureAtNodeCalc::declareTypes()
  {
    // Input: number density of species weighted by 2*pi/m
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: <H*f> weighted by 2*pi/m
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: <dH/dv*f> weighted by 2*pi/m
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Total temperature of species
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
