/**
 * @file	LcMaxwellEquation.cpp
 *
 * @brief	Maxwell equations for electromagnetism.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMaxwellEquation.h>

namespace Lucee
{
/** Class id: this is used by registration system */
  const char *MaxwellEquation::id = "Maxwell";

// constants to make indexing conserved variables easier
  static const unsigned EX = 0;
  static const unsigned EY = 1;
  static const unsigned EZ = 2;
  static const unsigned BX = 3;
  static const unsigned BY = 4;
  static const unsigned BZ = 5;

  MaxwellEquation::MaxwellEquation()
    : Lucee::HyperEquation(6, 3)
  { // 6 eqns and 3 waves: 0, c, -c
  }

  void
  MaxwellEquation::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::HyperEquation::readInput(tbl);
    if (tbl.hasNumber("lightSpeed"))
      lightSpeed = tbl.getNumber("lightSpeed");
    else
      throw Lucee::Except("MaxwellEquation::readInput: Must specify speed of light, 'lightSpeed'");
    lightSpeed2 = lightSpeed*lightSpeed;
  }

  void
  MaxwellEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    c.rotateVecToLocal(&inQ[0], &outQ[0]); // electric field
    c.rotateVecToLocal(&inQ[3], &outQ[3]); // magnetic field
  }

  void
  MaxwellEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    c.rotateVecToGlobal(&inQ[0], &outQ[0]); // electric field
    c.rotateVecToGlobal(&inQ[3], &outQ[3]); // magnetic field
  }

  void
  MaxwellEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    f[0] = 0.0;
    f[1] = lightSpeed2*q[BZ];
    f[2] = -lightSpeed2*q[BY];
    f[3] = 0.0;
    f[4] = -q[EZ];
    f[5] = q[EY];
  }

  void
  MaxwellEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    s[0] = -lightSpeed;
    s[1] = lightSpeed;
  }

  void
  MaxwellEquation::primitive(const double* q, double* v) const
  {
    for (unsigned i=0; i<6; ++i)
      v[i] = q[i];
  }

  void
  MaxwellEquation::conserved(const double* v, double* q) const
  {
    for (unsigned i=0; i<6; ++i)
      q[i] = v[i];
  }

  void
  MaxwellEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double c1 = 1/lightSpeed;
// project jump onto left eigenvectors (see Tech Note 1012)
    double a0 = jump[3];
    double a1 = jump[0];
    double a2 = 0.5*(c1*jump[1] + jump[5]);
    double a3 = 0.5*(-c1*jump[2] + jump[4]);
    double a4 = 0.5*(-c1*jump[1] + jump[5]);
    double a5 = 0.5*(c1*jump[2] + jump[4]);

// compute waves (see Tech Note 1012)

// wave 1: eigenvalue is 0 (multiplicity 2). We store the zero wave as
// this is needed when splitting jump in flux rather than in fields.
    waves(0,0) = a1;
    waves(1,0) = 0.0;
    waves(2,0) = 0.0;
    waves(3,0) = a0;
    waves(4,0) = 0.0;
    waves(5,0) = 0.0;
    s[0] = 0.0;

// wave 2: eigenvalue is c (multiplicity 2)
    waves(0,1) = 0.0;
    waves(1,1) = a2*lightSpeed;
    waves(2,1) = -a3*lightSpeed;
    waves(3,1) = 0.0;
    waves(4,1) = a3;
    waves(5,1) = a2;
    s[1] = lightSpeed;

// wave 2: eigenvalue is c (multiplicity 2)
    waves(0,2) = 0.0;
    waves(1,2) = -a4*lightSpeed;
    waves(2,2) = a5*lightSpeed;
    waves(3,2) = 0.0;
    waves(4,2) = a5;
    waves(5,2) = a4;
    s[2] = -lightSpeed;
  }
}
