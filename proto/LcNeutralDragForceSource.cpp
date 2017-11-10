/**
 * @file	LcNeutralDragForceSource.cpp
 *
 * @brief	Source for computing NeutralDrag force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcNeutralDragForceSource.h>

namespace Lucee
{
// set id for creators
  const char *NeutralDragForceSource::id = "NeutralDragForce";

  NeutralDragForceSource::NeutralDragForceSource()
    : Lucee::PointSourceIfc(4, 4)
  { 
    // takes in [rho, rho*u, rho*v, rho*w,] and computes
    // sources for [rho*u, rho*v, rho*w, Er]
  }

  void
  NeutralDragForceSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
    // get chargeIn and mass of species
    if (tbl.hasNumVec("velocityNeut"))
    {
      velocityNeut = tbl.getNumVec("velocityNeut");
      if (velocityNeut.size() != 3)
        throw Lucee::Except(
          "NeutralDragForce::readInput: 'velocityNeut' table must have size of 3");
    }
    else
      throw Lucee::Except("NeutralDragForce::readInput: 'velocityNeut' table must be specified");

    if (tbl.hasNumber("nu_ei"))
      nu_ei = tbl.getNumber("nu_ei");
    else
      throw Lucee::Except("NeutralDragForce::readInput: Electron-ion collision frequency 'nu_ei' must be specified");
    if (tbl.hasNumber("nu_en"))
      nu_en = tbl.getNumber("nu_en");
    else
      throw Lucee::Except("NeutralDragForce::readInput: Electron-neutral collision frequency 'nu_en' must be specified");
    if (tbl.hasNumber("nu_in"))
      nu_in = tbl.getNumber("nu_in");
    else
      throw Lucee::Except("NeutralDragForce::readInput: Ion-neutral collision frequency 'nu' must be specified");
  }

  void
  NeutralDragForceSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    double u_e, v_e, w_e;
    double u_i, v_i, w_i;

    double relU_ei, relV_ei, relW_ei;
    double relU_en, relV_en, relW_en;
    double relU_in, relV_in, relW_in;
    // takes in [rho_e, rho_e*u_e, rho_e*v_e, rho_e*w_e, e_e
    //           rho_i, rho_i*u_i, rho_i*v_i, rho_i*w_i, e_i]
    double rho_e = this->getData(0);
    double rhou_e = this->getData(1);
    double rhov_e = this->getData(2);
    double rhow_e = this->getData(3);
    double rho_i = this->getData(5);
    double rhou_i = this->getData(6);
    double rhov_i = this->getData(7);
    double rhow_i = this->getData(8);

    // get velocities
    u_e = rhou_e/rho_e;
    v_e = rhov_e/rho_e;
    w_e = rhow_e/rho_e;
    u_i = rhou_i/rho_i;
    v_i = rhov_i/rho_i;
    w_i = rhow_i/rho_i;

    // get raltive velocities
    relU_ei = u_e - u_i;
    relV_ei = v_e - v_i;
    relW_ei = w_e - w_i;
    relU_en = u_e - velocityNeut[0];
    relV_en = v_e - velocityNeut[1];
    relW_en = w_e - velocityNeut[2];
    relU_in = u_i - velocityNeut[0];
    relV_in = v_i - velocityNeut[1];
    relW_in = w_i - velocityNeut[2];

    // computes sources
    src[0] = 0;
    src[5] = 0;
    
    src[1] = - nu_ei*rho_e*relU_ei - nu_en*rho_e*relU_en;
    src[2] = - nu_ei*rho_e*relV_ei - nu_en*rho_e*relV_en;
    src[3] = - nu_ei*rho_e*relW_ei - nu_en*rho_e*relW_en;
    src[6] = nu_ei*rho_i*relU_ei - nu_in*rho_i*relU_in;
    src[7] = nu_ei*rho_i*relV_ei - nu_in*rho_i*relV_in;
    src[8] = nu_ei*rho_i*relW_ei - nu_in*rho_i*relW_in;
    
    src[4] = nu_ei*rho_e*(relU_ei*relU_ei + 
			  relV_ei*relV_ei + 
			  relW_ei*relW_ei) +
      nu_en*rho_e*(relU_en*relU_en + 
		   relV_en*relV_en + 
		   relW_en*relW_en);
    src[9] = -nu_ei*rho_i*(relU_ei*relU_ei + 
			   relV_ei*relV_ei + 
			   relW_ei*relW_ei) +
      nu_in*rho_i*(relU_in*relU_in + 
		   relV_in*relV_in + 
		   relW_in*relW_in);
  }

  void
  NeutralDragForceSource::getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac)
  {
  }
}
