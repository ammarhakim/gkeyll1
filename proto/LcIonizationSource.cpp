/**
 * @file	LcIonizationSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcIonizationSource.h>

namespace Lucee
{
  // set id for creators
  const char *IonizationSource::id = "Ionization";

  IonizationSource::IonizationSource()
    : Lucee::PointSourceIfc(13, 4)
  { 
    // takes in [rho, rho*u, rho*v, rho*w, Er, .., Bx, By, Bz] and computes
    // sources for [rho, Er]
    // for both electrons and ions
  }

  void
  IonizationSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
    // get ionization constant, density and gas gamma
    if (tbl.hasNumber("ionizationConst"))
      ionizationConst = tbl.getNumber("ionizationConst");
    else
      throw Lucee::Except("IonizationSource::readInput: Must specify ionization constant using 'ionizationConst'");

    if (tbl.hasNumber("massElc"))
      massElc = tbl.getNumber("massElc");
    else
      throw Lucee::Except("IonizationSource::readInput: Must specify electron mass using 'massElc'");
    if (tbl.hasNumber("massIon"))
      massIon = tbl.getNumber("massIon");
    else
      throw Lucee::Except("IonizationSource::readInput: Must specify ion mass using 'massIon'");

    ionizationConstDensityElc = ionizationConst;
    ionizationConstDensityIon = ionizationConst*massIon/massElc;
    ionizationConstEnergy = 0.5*ionizationConst/massElc;

    if (tbl.hasNumber("gasGamma"))
      gasGamma = tbl.getNumber("gasGamma");
    else
      throw Lucee::Except("IonizationSource::readInput: Must specify gas gamma using 'gasGamma'");

    if (tbl.hasNumber("permeability"))
    {
      double permeability = tbl.getNumber("permeability");
      invPermeability = 1/permeability;
    }
    else
      throw Lucee::Except("IonizationSource::readInput: Must specify permeability using 'permeability'");
  }

  void
  IonizationSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    double B2mu; // magnetic field magnitude squared divided by permitivity
    double vTerm2Elc; // electron thermal velocity
    double vTerm2Ion; // ion thermal velocity

    // inputs
    double rhoElc = this->getData(0);
    double rhouElc = this->getData(1);
    double rhovElc = this->getData(2);
    double rhowElc = this->getData(3);
    double ErElc = this->getData(4);
    double rhoIon = this->getData(5);
    double rhouIon = this->getData(6);
    double rhovIon = this->getData(7);
    double rhowIon = this->getData(8);
    double ErIon = this->getData(9);
    double Bx = this->getData(10);
    double By = this->getData(11);
    double Bz = this->getData(12);

    // calculated intermediate values
    B2mu = (Bx*Bx+By*By+Bz*Bz)*invPermeability;
    vTerm2Elc = (gasGamma-1)*(ErElc-0.5*((rhouElc*rhouElc+rhovElc*rhovElc+rhowElc*rhowElc)/rhoElc+B2mu))/rhoElc;
    vTerm2Ion = (gasGamma-1)*(ErIon-0.5*((rhouIon*rhouIon+rhovIon*rhovIon+rhowIon*rhowIon)/rhoIon+B2mu))/rhoIon;

    // computes sources
    src[0] = ionizationConstDensityElc*rhoElc; // electron density
    src[1] = ionizationConstEnergy*rhoElc*vTerm2Elc; // electron energy
    src[2] = ionizationConstDensityIon*rhoElc; // ion density
    src[3] = ionizationConstEnergy*rhoElc*vTerm2Ion; // ion energy
  }

  void
  IonizationSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
  }
}
