/**
 * @file	LcPhMaxwellEquation.cpp
 *
 * @brief	Perfectly-Hyperbolic Maxwell equations for electromagnetism.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPhMaxwellEquation.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *PhMaxwellEquation::id = "PhMaxwell";

// constants to make indexing conserved variables easier
  static const unsigned EX = 0;
  static const unsigned EY = 1;
  static const unsigned EZ = 2;
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;
  static const unsigned PHI = 6; // electric field correction
  static const unsigned PSI = 7; // magnetic field correction

  PhMaxwellEquation::PhMaxwellEquation()
    : Lucee::HyperEquation(8, 6)
  { // 8 eqns and 6 waves: -c*chi_m, c*chi_m, -c*chi_e, c*chi_e, c, c, -c, -c
  }

  void
  PhMaxwellEquation::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::HyperEquation::readInput(tbl);

    if (tbl.hasNumber("lightSpeed"))
      lightSpeed = tbl.getNumber("lightSpeed");
    else
      throw Lucee::Except("PhMaxwellEquation::readInput: Must specify speed of light, 'lightSpeed'");
    lightSpeed2 = lightSpeed*lightSpeed;

    if (tbl.hasNumber("elcErrorSpeedFactor"))
      chi_e = tbl.getNumber("elcErrorSpeedFactor");
    else
      chi_e = 0.0;

    if (tbl.hasNumber("mgnErrorSpeedFactor"))
      chi_m = tbl.getNumber("mgnErrorSpeedFactor");
    else
      chi_m = 0.0;
  }

  void
  PhMaxwellEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    c.rotateVecToLocal(&inQ[0], &outQ[0]); // electric field
    c.rotateVecToLocal(&inQ[3], &outQ[3]); // magnetic field
    outQ[6] = inQ[6]; // electric correction potential
    outQ[7] = inQ[7]; // magnetic correction potential
  }

  void
  PhMaxwellEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    c.rotateVecToGlobal(&inQ[0], &outQ[0]); // electric field
    c.rotateVecToGlobal(&inQ[3], &outQ[3]); // magnetic field
    outQ[6] = inQ[6]; // electric correction potential
    outQ[7] = inQ[7]; // magnetic correction potential
  }

  void
  PhMaxwellEquation::flux(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& f)
  {
    f[0] = chi_e*lightSpeed2*q[PHI];
    f[1] = lightSpeed2*q[BZ];
    f[2] = -lightSpeed2*q[BY];
    f[3] = chi_m*q[PHI];
    f[4] = -q[EZ];
    f[5] = q[EY];
    f[6] = chi_e*q[EX];
    f[7] = chi_m*lightSpeed2*q[BX];
  }

  void
  PhMaxwellEquation::speeds(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& q, double s[2])
  {
    if (chi_m<1.0 && chi_e<1.0)
    {
      s[0] = -lightSpeed;
      s[1] = lightSpeed;
    }
    else if (chi_m > chi_e)
    {
      s[0] = -chi_m*lightSpeed;
      s[1] = chi_m*lightSpeed;
    }
    else
    {
      s[0] = -chi_e*lightSpeed;
      s[1] = chi_e*lightSpeed;
    }
  }

  void
  PhMaxwellEquation::primitive(const Lucee::ConstFieldPtr<double>& q, Lucee::FieldPtr<double>& v) const
  {
    for (unsigned i=0; i<8; ++i)
      v[i] = q[i];
  }

  void
  PhMaxwellEquation::conserved(const Lucee::ConstFieldPtr<double>& v, Lucee::FieldPtr<double>& q) const
  {
    for (unsigned i=0; i<8; ++i)
      q[i] = v[i];
  }

  void
  PhMaxwellEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double c1 = 1/lightSpeed;
// project jump onto left eigenvectors (see Tech Note 1012)
    double a0 = 0.5*(jump[3] - c1*jump[7]);
    double a1 = 0.5*(jump[3] + c1*jump[7]);
    double a2 = 0.5*(jump[0] - lightSpeed*jump[6]);
    double a3 = 0.5*(jump[0] + lightSpeed*jump[6]);
    double a4 = 0.5*(c1*jump[1] + jump[5]);
    double a5 = 0.5*(-c1*jump[2] + jump[4]);
    double a6 = 0.5*(-c1*jump[1] + jump[5]);
    double a7 = 0.5*(c1*jump[2] + jump[4]);

// compute waves (see Tech Note 1012)

// wave 1: eigenvalue is -c*chi_m
    waves(0,0) = 0.0;
    waves(1,0) = 0.0;
    waves(2,0) = 0.0;
    waves(3,0) = a0;
    waves(4,0) = 0.0;
    waves(5,0) = 0.0;
    waves(6,0) = 0.0;
    waves(7,0) = -a0*lightSpeed;
    s[0] = -lightSpeed*chi_m;

// wave 2: eigenvalue is c*chi_m
    waves(0,1) = 0.0;
    waves(1,1) = 0.0;
    waves(2,1) = 0.0;
    waves(3,1) = a1;
    waves(4,1) = 0.0;
    waves(5,1) = 0.0;
    waves(6,1) = 0.0;
    waves(7,1) = a1*lightSpeed;
    s[1] = lightSpeed*chi_m;

// wave 3: eigenvalue is -c*chi_e
    waves(0,2) = a2;
    waves(1,2) = 0.0;
    waves(2,2) = 0.0;
    waves(3,2) = 0.0;
    waves(4,2) = 0.0;
    waves(5,2) = 0.0;
    waves(6,2) = -a2*c1;
    waves(7,2) = 0.0;
    s[2] = -lightSpeed*chi_e;

// wave 4: eigenvalue is c*chi_e
    waves(0,3) = a3;
    waves(1,3) = 0.0;
    waves(2,3) = 0.0;
    waves(3,3) = 0.0;
    waves(4,3) = 0.0;
    waves(5,3) = 0.0;
    waves(6,3) = a3*c1;
    waves(7,3) = 0.0;
    s[3] = lightSpeed*chi_e;

// wave 5: eigenvalue is c (multiplicity 2)
    waves(0,4) = 0.0;
    waves(1,4) = a4*lightSpeed;
    waves(2,4) = -a5*lightSpeed;
    waves(3,4) = 0.0;
    waves(4,4) = a5;
    waves(5,4) = a4;
    waves(6,4) = 0.0;
    waves(7,4) = 0.0;
    s[4] = lightSpeed;

// wave 6: eigenvalue is -c (multiplicity 2)
    waves(0,5) = 0.0;
    waves(1,5) = -a6*lightSpeed;
    waves(2,5) = a7*lightSpeed;
    waves(3,5) = 0.0;
    waves(4,5) = a7;
    waves(5,5) = a6;
    waves(6,5) = 0.0;
    waves(7,5) = 0.0;
    s[5] = -lightSpeed;
  }
}
