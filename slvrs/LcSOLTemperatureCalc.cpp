/**
 * @file	LcSOLTemperatureCalc.cpp
 *
 * @brief	Simple updater to compute temperature node-by-node for 5D SOL simulations
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSOLTemperatureCalc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  const char *SOLTemperatureCalc::id = "SOLTemperatureCalc";

  SOLTemperatureCalc::SOLTemperatureCalc()
  {
  }

  void
  SOLTemperatureCalc::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<3> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<3> >("basis");
    else
      throw Lucee::Except("SOLTemperatureCalc::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("speciesMass"))
      speciesMass = tbl.getNumber("speciesMass");
    else
      throw Lucee::Except("SOLTemperatureCalc::readInput: Must specify electron mass using 'speciesMass'");
  }

  void
  SOLTemperatureCalc::initialize()
  {
    UpdaterIfc::initialize();
  }

  Lucee::UpdaterStatus
  SOLTemperatureCalc::update(double t)
  {
    const Lucee::StructuredGridBase<3>& grid 
      = this->getGrid<Lucee::StructuredGridBase<3> >();

    // Input number density
    const Lucee::Field<3, double>& numDensIn = this->getInp<Lucee::Field<3, double> >(0);
    // Input first parallel velocity moment
    const Lucee::Field<3, double>& mom1dir3In = this->getInp<Lucee::Field<3, double> >(1);
    // Input second parallel velocity moment
    const Lucee::Field<3, double>& mom2dir3In = this->getInp<Lucee::Field<3, double> >(2);
    // Input first mu moment
    const Lucee::Field<3, double>& mom1dir4In = this->getInp<Lucee::Field<3, double> >(3);
    // Input magnetic field profile
    const Lucee::Field<3, double>& bFieldIn = this->getInp<Lucee::Field<3, double> >(4);
    // Output dynvector of electron temperatures at wall
    Lucee::Field<3, double>& temperatureOut = this->getOut<Lucee::Field<3, double> >(0);

    Lucee::ConstFieldPtr<double> numDensPtr = numDensIn.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1dir3Ptr = mom1dir3In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom2dir3Ptr = mom2dir3In.createConstPtr();
    Lucee::ConstFieldPtr<double> mom1dir4Ptr = mom1dir4In.createConstPtr();
    Lucee::ConstFieldPtr<double> bFieldPtr = bFieldIn.createConstPtr();
    
    Lucee::FieldPtr<double> temperaturePtr = temperatureOut.createPtr();

    unsigned nlocal = nodalBasis->getNumNodes();

    double cellCentroid[3];
    int idx[3];
    Lucee::Region<3, int> localRgn = grid.getLocalRegion();
    Lucee::RowMajorSequencer<3> seq(localRgn);

    // Loop over local region
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      numDensIn.setPtr(numDensPtr, idx);
      mom1dir3In.setPtr(mom1dir3Ptr, idx);
      mom2dir3In.setPtr(mom2dir3Ptr, idx);
      mom1dir4In.setPtr(mom1dir4Ptr, idx);
      bFieldIn.setPtr(bFieldPtr, idx);
      temperatureOut.setPtr(temperaturePtr, idx);

      for (int i = 0; i < nlocal; i++)
      {
        double tPara = speciesMass*(mom2dir3Ptr[i]/numDensPtr[i] - 
            mom1dir3Ptr[i]*mom1dir3Ptr[i]/(numDensPtr[i]*numDensPtr[i]));
        double tPerp = bFieldPtr[i]*mom1dir4Ptr[i]/numDensPtr[i];

        temperaturePtr[i] = (tPara + 2.0*tPerp)/3.0;
      }
    }

    // If there are three output fields, we need to 
    if (this->getNumOutVars() == 3)
    {
      Lucee::Field<3, double>& temperatureParaOut = this->getOut<Lucee::Field<3, double> >(1);
      Lucee::Field<3, double>& temperaturePerpOut = this->getOut<Lucee::Field<3, double> >(2);
      Lucee::FieldPtr<double> temperatureParaPtr = temperatureParaOut.createPtr();
      Lucee::FieldPtr<double> temperaturePerpPtr = temperaturePerpOut.createPtr();
      
      seq.reset();
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        numDensIn.setPtr(numDensPtr, idx);
        mom1dir3In.setPtr(mom1dir3Ptr, idx);
        mom2dir3In.setPtr(mom2dir3Ptr, idx);
        mom1dir4In.setPtr(mom1dir4Ptr, idx);
        bFieldIn.setPtr(bFieldPtr, idx);

        temperatureParaOut.setPtr(temperatureParaPtr, idx);
        temperaturePerpOut.setPtr(temperaturePerpPtr, idx);

        for (int i = 0; i < nlocal; i++)
        {
          double tPara = speciesMass*(mom2dir3Ptr[i]/numDensPtr[i] - 
              mom1dir3Ptr[i]*mom1dir3Ptr[i]/(numDensPtr[i]*numDensPtr[i]));
          double tPerp = bFieldPtr[i]*mom1dir4Ptr[i]/numDensPtr[i];

          temperatureParaPtr[i] = tPara;
          temperaturePerpPtr[i] = tPerp;
        }
      }
    }
   

    return Lucee::UpdaterStatus();
  }

  void
  SOLTemperatureCalc::declareTypes()
  {
    // Input: number density of species
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: first parallel velocity moment of species
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: second parallel velocity moment of species
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: first mu moment of species
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Input: magnetic field profile
    this->appendInpVarType(typeid(Lucee::Field<3, double>));
    // Output: Total temperature of species
    this->appendOutVarType(typeid(Lucee::Field<3, double>));
  }
}
