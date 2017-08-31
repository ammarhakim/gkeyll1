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
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify ionization energy using 'phi'");
    if (tbl.hasNumber("A"))
      voronovA = tbl.getNumber("A");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify Voronov ionization fitting parameter A using 'A'");
    if (tbl.hasNumber("P"))
      voronovP = tbl.getNumber("P");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify Voronov ionization fitting parameter P using 'P'");
    if (tbl.hasNumber("X"))
      voronovX = tbl.getNumber("X");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify Voronov ionization fitting parameter X using 'X'");
    if (tbl.hasNumber("K"))
      voronovK = tbl.getNumber("K");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify Voronov ionization fitting parameter K using 'K'");

    if (tbl.hasNumber("massElc"))
      mass_e = tbl.getNumber("massElc");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify electron mass using 'massElc'");
    if (tbl.hasNumber("massIon"))
      mass_i = tbl.getNumber("massIon");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify ion mass using 'massIon'");
    if (tbl.hasNumber("massNeut"))
      mass_n = tbl.getNumber("massNeut");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify neutral mass using 'massNeut'");

    if (tbl.hasNumber("chargeElc"))
      charge_e = tbl.getNumber("chargeElc");
    else
      throw Lucee::Except("ThreeFluidCollisionSource::readInput: Must specify electron charge using 'chargeElc'");
    if (tbl.hasNumber("chargeIon"))
      charge_i = tbl.getNumber("chargeIon");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify ion charge using 'chargeIon'");

    if (tbl.hasNumber("gasGamma"))
      gasGamma = tbl.getNumber("gasGamma");
    else
      throw Lucee::Except("ThreeFluidCollision::readInput: Must specify gas gamma using 'gasGamma'");
  }

  void
  ThreeFluidCollisionSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    // inputs
    double rho_e = this->getData(0);
    double rhoUx_e = this->getData(1);
    double rhoUy_e = this->getData(2);
    double rhoUz_e = this->getData(3);
    double Edens_e = this->getData(4);
    double rho_i = this->getData(5);
    double rhoUx_i = this->getData(6);
    double rhoUy_i = this->getData(7);
    double rhoUz_i = this->getData(8);
    double Edens_i = this->getData(9);
    double rho_n = this->getData(10);
    double rhoUx_n = this->getData(11);
    double rhoUy_n = this->getData(12);
    double rhoUz_n = this->getData(13);
    double Edens_n = this->getData(14);

    double elemCharge, Z;
    double n_e, n_i, n_n;
    double T_e, T_i, T_n;
    double EdensK_e, EdensK_i, EdensK_n, EK_n;
    double ux_n, uy_n, uz_n;

    double voronovU, ionizationRate, gammaIonization;

    double tau_e;
    double sigma;
    double jx, jy, jz;
    double coulombLog;
    double Rx_ei, Ry_ei, Rz_ei;

    double Qion_n, Q_ie, Q_ei;

    elemCharge = fabs(charge_e);
    Z = charge_i / elemCharge;

    n_e = rho_e / mass_e;
    n_i = rho_i / mass_i;
    n_n = rho_n / mass_n;

    ux_n = rhoUx_n / rho_n;
    uy_n = rhoUy_n / rho_n;
    uz_n = rhoUz_n / rho_n;

    EdensK_e = 0.5*(rhoUx_e*rhoUx_e + rhoUy_e*rhoUy_e +
		    rhoUz_e*rhoUz_e)/rho_e;
    EdensK_i = 0.5*(rhoUx_i*rhoUx_i + rhoUy_i*rhoUy_i +
		    rhoUz_i*rhoUz_i)/rho_i;
    EdensK_n = 0.5*(rhoUx_n*rhoUx_n + rhoUy_n*rhoUy_n +
		    rhoUz_n*rhoUz_n)/rho_n;
    EK_n = EdensK_n*mass_n/rho_n;
    
    // all the following formulas are for temperatures in eV but
    // sumulations are running in SI; therefore elemCharge can
    // unexpectedly appear in the formulas
    T_e = (gasGamma-1)*(Edens_e - EdensK_e)/n_e / elemCharge;
    T_i = (gasGamma-1)*(Edens_i - EdensK_i)/n_i / elemCharge;
    T_n = (gasGamma-1)*(Edens_n - EdensK_n)/n_n / elemCharge;

    voronovU = phi / T_e;
    ionizationRate = voronovA*1e-6* (1 + voronovP*sqrt(voronovU)) / (voronovX + voronovU) * pow(voronovU, voronovK) * exp(-1.0*voronovU);
    gammaIonization = n_e*n_n*ionizationRate;

    if (T_e < 50) {
      coulombLog = 23.4 - 1.15*log(n_e*1e-6) + 3.45*log(T_e);
    }
    else {
      coulombLog = 25.3 - 1.15*log(n_e*1e-6) + 2.3*log(T_e);
    }

    // From Braginskii:
    tau_e = 3.5e11/coulombLog * sqrt(T_e*T_e*T_e)/Z * mass_e/rho_e;
    sigma = 1.96 * n_e * elemCharge*elemCharge * tau_e / mass_e;
    jx = charge_e*rhoUx_e/mass_e + charge_i*rhoUx_i/mass_i;
    jy = charge_e*rhoUy_e/mass_e + charge_i*rhoUy_i/mass_i;
    jz = charge_e*rhoUz_e/mass_e + charge_i*rhoUz_i/mass_i;
    Rx_ei = elemCharge*n_e*jx/sigma;
    Ry_ei = elemCharge*n_e*jy/sigma;
    Rz_ei = elemCharge*n_e*jz/sigma;

    Qion_n = 1.5*gammaIonization*T_n * elemCharge;
    Q_ie = 3*rho_e/mass_i/tau_e*(T_e-T_i) * elemCharge;
    Q_ei = (jx*jx + jy*jy + jz*jz)/sigma * elemCharge - Q_ie;
      
    // continuity sources
    src[0] = mass_e*gammaIonization;
    src[5] = mass_i*gammaIonization;
    src[10] = -mass_n*gammaIonization;

    // momentum sources
    src[1] = gammaIonization*mass_e*ux_n + Rx_ei;
    src[2] = gammaIonization*mass_e*uy_n + Ry_ei;
    src[3] = gammaIonization*mass_e*uz_n + Rz_ei;
    src[6] = gammaIonization*mass_i*ux_n;// - Rx_ei;
    src[7] = gammaIonization*mass_i*uy_n - Ry_ei;
    src[8] = gammaIonization*mass_i*uz_n - Rz_ei;
    src[11] = -gammaIonization*mass_n*ux_n; 
    src[12] = -gammaIonization*mass_n*uy_n;
    src[13] = -gammaIonization*mass_n*uz_n;

    // energy sources
    src[4] = (Rx_ei*rhoUx_e + Ry_ei*rhoUy_e + Rz_ei*rhoUz_e)/rho_e + 
      mass_e/mass_n*(gammaIonization*EK_n + Qion_n) + Q_ei; 
    src[9] = -(Rx_ei*rhoUx_i + Ry_ei*rhoUy_i + Rz_ei*rhoUz_i)/rho_i + 
      mass_i/mass_n*(gammaIonization*EK_n + Qion_n) + Q_ie;
    src[14] = -(gammaIonization*EK_n + Qion_n); 
  }

  void
  ThreeFluidCollisionSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
  }
}
