/**
 * @file	LcThreeFluidCollisionSource.cpp
 *
 * @brief	Source for three fluid collisions [Meier & Shumlak, 2012]
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcThreeFluidCollisionSource.h>

namespace Lucee
{
  // set id for creators
  const char *ThreeFluidCollisionSource::id = "ThreeFluidCollision";

  ThreeFluidCollisionSource::ThreeFluidCollisionSource()
    : Lucee::PointSourceIfc(15, 15)
  { 
    // takes in [rho, rho*u, rho*v, rho*w, Er] for electrons, ions, 
    // and neutrals and computes the corresponding sources 
  }

  void
  ThreeFluidCollisionSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
    // get Voronov fit parameters
    if (tbl.hasNumber("phi"))
      phi = tbl.getNumber("phi");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify ionization energy using 'phi'");
    if (tbl.hasNumber("A"))
      voronovA = tbl.getNumber("A");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify Voronov ionization fitting parameter A using 'A'");
    if (tbl.hasNumber("P"))
      voronovP = tbl.getNumber("P");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify Voronov ionization fitting parameter P using 'P'");
    if (tbl.hasNumber("X"))
      voronovX = tbl.getNumber("X");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify Voronov ionization fitting parameter X using 'X'");
    if (tbl.hasNumber("K"))
      voronovK = tbl.getNumber("K");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify Voronov ionization fitting parameter K using 'K'");

    if (tbl.hasNumber("massElc"))
      massElc = tbl.getNumber("massElc");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify electron mass using 'massElc'");
    if (tbl.hasNumber("massIon"))
      massIon = tbl.getNumber("massIon");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify ion mass using 'massIon'");
    if (tbl.hasNumber("massNeut"))
      massIon = tbl.getNumber("massNeut");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify neutrals mass using 'massNeut'");

    if (tbl.hasNumber("chargeElc"))
      chargeElc = tbl.getNumber("chargeElc");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify electron charge using 'chargeElc'");
    if (tbl.hasNumber("chargeIon"))
      chargeIon = tbl.getNumber("chargeIon");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify ion charge using 'chargeIon'");


    if (tbl.hasNumber("gasGamma"))
      gasGamma = tbl.getNumber("gasGamma");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify gas gamma using 'gasGamma'");
  }

  void
  ThreeFluidCollisionSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    // inputs
    double rhoElc = this->getData(0);
    double rhoUxElc = this->getData(1);
    double rhoUyElc = this->getData(2);
    double rhoUzElc = this->getData(3);
    double ErElc = this->getData(4);
    double rhoIon = this->getData(5);
    double rhoUxIon = this->getData(6);
    double rhoUyIon = this->getData(7);
    double rhoUzIon = this->getData(8);
    double ErIon = this->getData(9);
    double rhoNeut = this->getData(10);
    double rhoUxNeut = this->getData(11);
    double rhoUyNeut = this->getData(12);
    double rhoUzNeut = this->getData(13);
    double ErNeut = this->getData(14);

    double elemCharge, Z;
    double invRhoElc, invRhoIon, invRhoNeut;
    double invMassElc, invMassIon, invMassNeut;
    double nElc, nNeut;
    double TElc, TIon, TNeut;
    double EkinElc, EkinIon, EkinNeut;
    double uxNeut, uyNeut, uzNeut;

    double voronovU, ionizationRate, gammaIonization;

    double tauElc;
    double sigma;
    double jx, jy, jz;
    double coulombLog;
    double Reix, Reiy, Reiz;

    double QnIonization, Qie, Qei;

    elemCharge = abs(chargeElc);
    Z = chargeIon/elemCharge;
    invRhoElc = 1.0/rhoElc;
    invRhoIon = 1.0/rhoIon;
    invRhoNeut = 1.0/rhoNeut;
    invMassElc = 1.0/massElc;
    invMassIon = 1.0/massIon;
    invMassNeut = 1.0/massNeut;

    nElc = rhoElc*invMassElc;
    nNeut = rhoNeut*invMassNeut;

    EkinElc = 0.5*invRhoElc*(rhoUxElc*rhoUxElc + rhoUyElc*rhoUyElc +
			     rhoUzElc*rhoUzElc);
    EkinIon = 0.5*invRhoIon*(rhoUxIon*rhoUxIon + rhoUyIon*rhoUyIon +
			     rhoUzIon*rhoUzIon);
    EkinNeut = 0.5*invRhoNeut*(rhoUxNeut*rhoUxNeut + rhoUyNeut*rhoUyNeut +
			       rhoUzNeut*rhoUzNeut);
    
    TElc = (gasGamma-1)*(ErElc - EkinElc)*invRhoElc / elemCharge; // in eV
    TIon = (gasGamma-1)*(ErIon - EkinIon)*invRhoIon / elemCharge; // in eV
    TNeut = (gasGamma-1)*(ErNeut - EkinNeut)*invRhoNeut / elemCharge;
    uxNeut = rhoUxNeut*invMassNeut;
    uyNeut = rhoUyNeut*invMassNeut;
    uzNeut = rhoUzNeut*invMassNeut;
    voronovU = phi/TElc;
    ionizationRate = voronovA*1e-6*(1 + voronovP*sqrt(voronovU))/(voronovX + voronovU)*pow(voronovU, voronovK)*exp(-1.0*voronovU);
    gammaIonization = nElc*nNeut*ionizationRate;

    if (TElc < 50) {
      coulombLog = 23.4 - 1.15*log(nElc/1e6) + 3.45*log(TElc);
    }
    else {
      coulombLog = 25.3 - 1.15*log(nElc/1e6) + 2.3*log(TElc);
    }

    tauElc = 3.5e11/coulombLog*sqrt(TElc*TElc*TElc)/Z * massElc*invRhoElc;
    sigma = 1.96*nElc*elemCharge*elemCharge*tauElc*invMassElc;
    jx = chargeElc*rhoUxElc*invMassElc + chargeIon*rhoUxIon*invMassIon;
    jy = chargeElc*rhoUyElc*invMassElc + chargeIon*rhoUyIon*invMassIon;
    jz = chargeElc*rhoUzElc*invMassElc + chargeIon*rhoUzIon*invMassIon;
    Reix = elemCharge*nElc*jx/sigma;
    Reiy = elemCharge*nElc*jy/sigma;
    Reiz = elemCharge*nElc*jz/sigma;

    QnIonization = 1.5*gammaIonization*TNeut;
    Qie = 3*rhoElc*invMassIon/tauElc*(TElc-TIon);
    Qei = (jx*jx + jy*jy + jz*jz)/sigma - Qie;
      
    // continuity sources
    src[0] = massElc*gammaIonization;
    src[5] = massIon*gammaIonization;
    src[10] = -1.0*massNeut*gammaIonization;

    // momentum sources
    src[1] = gammaIonization*massElc*uxNeut + Reix;
    src[2] = gammaIonization*massElc*uyNeut + Reiy;
    src[3] = gammaIonization*massElc*uzNeut + Reiz;
    src[6] = gammaIonization*massIon*uxNeut - Reix;
    src[7] = gammaIonization*massIon*uyNeut - Reiy;
    src[8] = gammaIonization*massIon*uzNeut - Reiz;
    src[11] = -1.0*gammaIonization*rhoUxNeut;
    src[12] = -1.0*gammaIonization*rhoUyNeut;
    src[13] = -1.0*gammaIonization*rhoUzNeut;

    // energy sources
    src[4] = massElc*invMassNeut*(gammaIonization*EkinNeut + QnIonization)
      + Qei; 
    src[9] = massIon*invMassNeut*(gammaIonization*EkinNeut + QnIonization)
      + Qie;
    src[14] = -1.0*(gammaIonization*EkinNeut + QnIonization); 
  }

  void
  ThreeFluidCollisionSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
  }
}
