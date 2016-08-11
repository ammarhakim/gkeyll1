/**
 * @file	LcIsoThermEulerEquation.cpp
 *
 * @brief	Euler equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcIsoThermEulerEquation.h>
#include <LcMathLib.h>

// std includes
#include <cmath>

namespace Lucee
{
// set ids for creators
  const char *IsoThermEulerEquation::id = "IsothermalEuler";

  IsoThermEulerEquation::IsoThermEulerEquation()
    : Lucee::HyperEquation(4, 3)
  {
  }
  
  void
  IsoThermEulerEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);

    if (tbl.hasNumber("thermalVelocity"))
      vt = tbl.getNumber("thermalVelocity");
    else
      throw Lucee::Except("IsoThermEulerEquation::readInput: Must specify thermalVelocity");

    correct = false; // by default do not correct
    if (tbl.hasBool("correct"))
      correct = tbl.getBool("correct");

    minDensity = 0.0;
    if (tbl.hasNumber("minDensity"))
      minDensity = tbl.getNumber("minDensity");

    numFlux = NF_ROE;
    if (tbl.hasString("numericalFlux"))
    {
      std::string nf = tbl.getString("numericalFlux");
      if (nf == "roe")
        numFlux = NF_ROE;
      else if (nf == "lax")
        numFlux = NF_LAX;
      else
      {
        Lucee::Except lce("IsoThermEulerEquation::readInput: 'numericalFlux' ");
        lce << nf << " not recognized!" << std::endl;
          throw lce;
      }
    }

    if (numFlux == NF_LAX)
    {
      useIntermediateWave = false;
      if (tbl.hasBool("useIntermediateWave"))
// this flag activates the intermediate state wave: it is best to use
// this when using Lax fluxes with second order wave-propagation
// scheme so as to avoid asymmetries.
        useIntermediateWave = tbl.getBool("useIntermediateWave");

// adjust number of waves accordingly
      if (useIntermediateWave)
        this->setNumWaves(2);
      else
        this->setNumWaves(1);
    }
  }

  void
  IsoThermEulerEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
  }

  void
  IsoThermEulerEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
  }

  void
  IsoThermEulerEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double rho = getSafeRho(q[0]);
    double rho1 = 1/rho;

    f[0] = q[1];
    f[1] = rho1*q[1]*q[1] + vt*vt*rho;
    f[2] = rho1*q[1]*q[2];
    f[3] = rho1*q[1]*q[3];
  }

  void
  IsoThermEulerEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    double rho = getSafeRho(q[0]);
    double u = q[1]/rho;
    s[0] = u-vt;
    s[1] = u+vt;
  }

  double
  IsoThermEulerEquation::maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q)
  {
    double rho = getSafeRho(q[0]);
    double u = q[1]/rho;
    return std::fabs(u)+vt;
  }

  void
  IsoThermEulerEquation::primitive(const double* q, double* v) const
  {
    double rho = getSafeRho(q[0]);
    v[0] = rho;
    v[1] = q[1]/rho;
    v[2] = q[2]/rho;
    v[3] = q[3]/rho;
  }

  void
  IsoThermEulerEquation::conserved(const double* v, double* q) const
  {
    q[0] = v[0];
    q[1] = v[0]*v[1];
    q[2] = v[0]*v[2];
    q[3] = v[0]*v[3];
  }

  void
  IsoThermEulerEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,    
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    if (numFlux == NF_ROE)
      wavesRoe(c, jump, ql, qr, waves, s);
    else if (numFlux == NF_LAX)
      wavesLax(c, jump, ql, qr, waves, s);
    else
    { /* this can't happen */ }
  }

  void
  IsoThermEulerEquation::wavesRoe(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double rhol = getSafeRho(ql[0]), rhor = getSafeRho(qr[0]);
// compute Roe averages for use in Riemann solver
    double rhsqrtl = std::sqrt(rhol); // left sqrt(rho)
    double rhsqrtr = std::sqrt(rhor); // right sqrt(rho)

    double rhsq2 = rhsqrtl + rhsqrtr;
// compute Roe averaged velocity components
    double u = (ql[1]/rhsqrtl  + qr[1]/rhsqrtr)/rhsq2;
    double v = (ql[2]/rhsqrtl  + qr[2]/rhsqrtr)/rhsq2;
    double w = (ql[3]/rhsqrtl  + qr[3]/rhsqrtr)/rhsq2;

// project onto left eigenvectors with Roe averaged values (see isothermal-euler.mac)
    double a0 = jump[0]*(vt+u)/vt/2.0-jump[1]/vt/2.0;
    double a1 = jump[2]-jump[0]*v;
    double a2 = jump[3]-jump[0]*w;
    double a3 = jump[0]*(vt-u)/vt/2.0+jump[1]/vt/2.0;

// compute waves

// wave 1: eigenvalue is u-c
    waves(0,0) = a0;
    waves(1,0) = a0*(u-vt);
    waves(2,0) = a0*v;
    waves(3,0) = a0*w;
    s[0] = u-vt;

// wave 2: 2 eigenvectors for the repeated roots are combined into a
// single wave
    waves(0,1) = 0;
    waves(1,1) = 0;
    waves(2,1) = a1;
    waves(3,1) = a2;
    s[1] = u;

// wave 3: eigenvalue is u+c
    waves(0,2) = a3;
    waves(1,2) = a3*(u+vt);
    waves(2,2) = a3*v;
    waves(3,2) = a3*w;
    s[2] = u+vt;
  }

  void
  IsoThermEulerEquation::wavesLax(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    if (useIntermediateWave)
    {
// this uses the HLLE technique to construct an intermediate state
      std::vector<const double*> auxVars;
      double fl[4], fr[4];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[0]+sr[0]);
      s[1] = 0.5*(sl[1]+sr[1]);

// compute intermediate HLLE state
      double sdiff1 = 1/(s[1]-s[0]);
      double qHHLE[4];
      for (unsigned i=0; i<4; ++i)
        qHHLE[i] = (s[1]*qr[i]-s[0]*ql[i]+fl[i]-fr[i])*sdiff1;

// compute waves
      for (unsigned i=0; i<4; ++i)
      {
        waves(i,0) = qHHLE[i]-ql[i];
        waves(i,1) = qr[i]-qHHLE[i];
      }
    }
    else
    {
// use a single wave
      for (unsigned i=0; i<4; ++i)
        waves(i,0) = qr[i]-ql[i];
      
      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[1]+sr[1]);
    }
  }

  double
  IsoThermEulerEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
// NOTE: This numerical flux is using Lax-Fluxes

    double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));

    double fl[4], fr[4];
    flux(c, ql, auxVarsl, fl);
    flux(c, qr, auxVarsr, fr);

    for (unsigned i=0; i<4; ++i)
      f[i] = 0.5*(fr[i]+fl[i]) - 0.5*absMaxs*(qr[i]-ql[i]);

    return absMaxs;
  }

  void
  IsoThermEulerEquation::qFluctuations(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s,
    Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq)
  {
    if (numFlux == NF_ROE)
      Lucee::HyperEquation::qFluctuations(c, ql, qr, waves, s, amdq, apdq);
    else if (numFlux == NF_LAX)
    {
// This might appear strange: we are arranging the fluctuations to
// effectively give Lax fluxes. Note that the sum of the fluctuations
// must be the jump in the physical flux, in all cases, as it does
// below.
      std::vector<const double*> auxVars;
      double fl[4], fr[4];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));
      double amdqL[4], apdqL[4];
      for (unsigned m=0; m<4; ++m)
      {
        amdqL[m] = 0.5*(fr[m]-fl[m] - absMaxs*(qr[m]-ql[m]));
        apdqL[m] = 0.5*(fr[m]-fl[m] + absMaxs*(qr[m]-ql[m]));
      }

// These rotations to global coordinate system are needed as the
// solvers expect fluctuations in global coordinates. This is contrary
// to most other methods in this class which work in the local
// coordinate system.
      rotateToGlobal(c, amdqL, &amdq[0]);
      rotateToGlobal(c, apdqL, &apdq[0]);
    }
  }

  bool
  IsoThermEulerEquation::isInvariantDomain(const double* q) const
  {
    double rho = q[0];
    if (Lucee::isNan(rho) || (rho<=0.0))
      return false;
    return true;
  }
      
  double
  IsoThermEulerEquation::getSafeRho(double rho) const
  {
    if (correct && (rho<minDensity))
      return minDensity;
    return rho;
  }
}
