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
#include <LcMathLib.h>
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
  { // 8 eqns and 6 waves: -c*chi_m, c*chi_m, -c*chi_e, c*chi_e, -c, -c, c, c
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

// store maximum speed for future use
    maxWaveSpeed = Lucee::max3(lightSpeed, chi_e*lightSpeed, chi_m*lightSpeed);

    numFlux = NF_UPWIND;
    if (tbl.hasString("numericalFlux"))
    {
      std::string nf = tbl.getString("numericalFlux");
      if (nf == "central")
        numFlux = NF_CENTRAL;
      else if (nf == "lax")
        numFlux = NF_LAX;
      else if (nf == "upwind")
        numFlux = NF_UPWIND;
      else
      {
        Lucee::Except lce("PhMaxwellEquation::readInput: 'numericalFlux' ");
        lce << nf << " not recognized!" << std::endl;
          throw lce;
      }
    }

// flag to indicate if there is a static field
    hasStatic = false;
    if (tbl.hasBool("hasStaticField"))
      hasStatic = tbl.getBool("hasStaticField");
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
  PhMaxwellEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
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
  PhMaxwellEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
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
  PhMaxwellEquation::primitive(const double* q, double* v) const
  {
    for (unsigned i=0; i<8; ++i)
      v[i] = q[i];
  }

  void
  PhMaxwellEquation::conserved(const double* v, double* q) const
  {
    for (unsigned i=0; i<8; ++i)
      q[i] = v[i];
  }

  void
  PhMaxwellEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,    
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double c1 = 1/lightSpeed;

// total jump may include static EM fields passed as aux variable    
    double totalJump[8];
    for (unsigned i=0; i<8; ++i)
      totalJump[i] = jump[i];
// add in extra jump if we have static fields    
    if (hasStatic)
    {
      double emL[8], emR[8];
      this->rotateToLocal(c, auxVarsl[0], emL);
      this->rotateToLocal(c, auxVarsr[0], emR);
      totalJump[BX] += emR[BX]-emL[BX];
    }
    
// project jump onto left eigenvectors (see Tech Note 1012)
    double a0 = 0.5*(totalJump[3] - c1*totalJump[7]);
    double a1 = 0.5*(totalJump[3] + c1*totalJump[7]);
    double a2 = 0.5*(totalJump[0] - lightSpeed*totalJump[6]);
    double a3 = 0.5*(totalJump[0] + lightSpeed*totalJump[6]);
    double a4 = 0.5*(totalJump[1] - lightSpeed*totalJump[5]);
    double a5 = 0.5*(totalJump[2] + lightSpeed*totalJump[4]);
    double a6 = 0.5*(totalJump[1] + lightSpeed*totalJump[5]);
    double a7 = 0.5*(totalJump[2] - lightSpeed*totalJump[4]);

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

// wave 5: eigenvalue is -c (multiplicity 2)
    waves(0,4) = 0.0;
    waves(1,4) = a4;
    waves(2,4) = a5;
    waves(3,4) = 0.0;
    waves(4,4) = a5*c1;
    waves(5,4) = -a4*c1;
    waves(6,4) = 0.0;
    waves(7,4) = 0.0;
    s[4] = -lightSpeed;

// wave 6: eigenvalue is c (multiplicity 2)
    waves(0,5) = 0.0;
    waves(1,5) = a6;
    waves(2,5) = a7;
    waves(3,5) = 0.0;
    waves(4,5) = -a7*c1;
    waves(5,5) = a6*c1;
    waves(6,5) = 0.0;
    waves(7,5) = 0.0;
    s[5] = lightSpeed;
  }

  double
  PhMaxwellEquation::maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q)
  {
    return maxWaveSpeed;
  }

  double 
  PhMaxwellEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
#define AVG(a,b) 0.5*(a+b)
#define JMP(a,b) 0.5*(a-b)

    double absMaxs = maxWaveSpeed;
    double qu[8], fl[8], fr[8];

    if (numFlux == NF_CENTRAL)
    {
      flux(c, ql, auxVarsl, fl);
      flux(c, qr, auxVarsr, fr);
      for (unsigned i=0; i<8; ++i)
        f[i] = 0.5*(fr[i]+fl[i]);
    }
    else if (numFlux == NF_LAX)
    {
      flux(c, ql, auxVarsl, fl);
      flux(c, qr, auxVarsr, fr);
      for (unsigned i=0; i<8; ++i)
        f[i] = 0.5*(fr[i]+fl[i]) - 0.5*absMaxs*(qr[i]-ql[i]);
    }
    else if (numFlux == NF_UPWIND)
    {
// first compute upwinded field at interface
      qu[EX] = AVG(qr[EX],ql[EX]) - lightSpeed*JMP(qr[PHI],ql[PHI]);
      qu[EY] = AVG(qr[EY],ql[EY]) - lightSpeed*JMP(qr[BZ],ql[BZ]);
      qu[EZ] = AVG(qr[EZ],ql[EZ]) + lightSpeed*JMP(qr[BY],ql[BY]);

      qu[BX] = AVG(qr[BX],ql[BX]) - JMP(qr[PSI],ql[PSI])/lightSpeed;
      qu[BY] = AVG(qr[BY],ql[BY]) + JMP(qr[EZ],ql[EZ])/lightSpeed;
      qu[BZ] = AVG(qr[BZ],ql[BZ]) - JMP(qr[EY],ql[EY])/lightSpeed;

      qu[PHI] = AVG(qr[PHI],ql[PHI]) - JMP(qr[EX],ql[EX])/lightSpeed;
      qu[PSI] = AVG(qr[PSI],ql[PSI]) - lightSpeed*JMP(qr[BX],ql[BX]);

// now compute interface flux: this is numerical flux
      flux(c, qu, auxVarsr, f);
    }
    else { /* Can't happen */ }

    return absMaxs;
#undef AVG
#undef JMP
  }

  bool
  PhMaxwellEquation::isInvariantDomain(const double* q) const
  {
    return true;
  }
}
