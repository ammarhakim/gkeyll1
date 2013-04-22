/**
 * @file	LcEulerEquation.cpp
 *
 * @brief	Euler equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEulerEquation.h>

// std includes
#include <cmath>

namespace Lucee
{
// set ids for creators
  const char *EulerEquation::id = "Euler";

  EulerEquation::EulerEquation()
    : Lucee::HyperEquation(5, 3)
  {
  }
  
  void
  EulerEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);
    if (tbl.hasNumber("gasGamma"))
      gas_gamma = tbl.getNumber("gasGamma");
    else
      throw Lucee::Except("EulerEquation::readInput: Must specify gasGamma (gas constant)");

    correct = false; // by default do not correct
    if (tbl.hasBool("correct"))
      correct = tbl.getBool("correct");

    minDensity = 0.0;
    if (tbl.hasNumber("minDensity"))
      minDensity = tbl.getNumber("minDensity");

    minPressure = 0.0;
    if (tbl.hasNumber("minPressure"))
      minPressure = tbl.getNumber("minPressure");

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
        Lucee::Except lce("EulerEquation::readInput: 'numericalFlux' ");
        lce << nf << " not recognized!" << std::endl;
          throw lce;
      }
    }

    if (numFlux == NF_LAX)
// for Lax fluxes there is just a single wave
      this->setNumWaves(1);
  }

  void
  EulerEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
    outQ[4] = inQ[4]; // energy density
  }

  void
  EulerEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
    outQ[4] = inQ[4]; // energy density
  }

  void
  EulerEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    double rho1 = 1/getSafeRho(q[0]);
    double pr = pressure(q);

    f[0] = q[1]; // rho*u
    f[1] = rho1*q[1]*q[1] + pr; // rho*u*u + pr
    f[2] = rho1*q[1]*q[2]; // rho*u*v
    f[3] = rho1*q[1]*q[3]; // rho*u*w
    f[4] = rho1*(q[4]+pr)*q[1]; // (E+pr)*u
  }

  void
  EulerEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    double rho = getSafeRho(q[0]);
    double pr = pressure(q);
    double cs = std::sqrt(gas_gamma*pr/rho); // sound speed
    double u = q[1]/rho; // fluid velocity
    s[0] = u-cs;
    s[1] = u+cs;
  }

  double
  EulerEquation::maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q)
  {
    double rho1 = 1/getSafeRho(q[0]);
    double pr = pressure(&q[0]);
    double cs = std::sqrt(gas_gamma*pr*rho1); // sound speed
    double u = q[1]*rho1; // fluid velocity
    return std::fabs(u)+cs;
  }

  void
  EulerEquation::primitive(const double* q, double* v) const
  {
    double rho = getSafeRho(q[0]);
    v[0] = rho; // rho
    v[1] = q[1]/rho; // u
    v[2] = q[2]/rho; // v
    v[3] = q[3]/rho; // w
    v[4] = pressure(q); // pressure
  }

  void
  EulerEquation::conserved(const double* v, double* q) const
  {
    q[0] = v[0]; // rho
    q[1] = v[0]*v[1]; // rho*u
    q[2] = v[0]*v[2]; // rho*v
    q[3] = v[0]*v[3]; // rho*w
    q[4] = v[4]/(gas_gamma-1) + 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]); // energy
  }

  void
  EulerEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
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
  EulerEquation::wavesRoe(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    double gas_gamma1 = gas_gamma-1;

    double rhol = getSafeRho(ql[0]), rhor = getSafeRho(qr[0]);
// compute Roe averages for use in Riemann solver
    double rhsqrtl = std::sqrt(rhol); // left sqrt(rho)
    double rhsqrtr = std::sqrt(rhor); // right sqrt(rho)
    double pl = pressure(&ql[0]); // left pressure
    double pr = pressure(&qr[0]); // right pressure

    double rhsq2 = rhsqrtl + rhsqrtr;
// compute Roe averaged velocity components
    double u = (ql[1]/rhsqrtl  + qr[1]/rhsqrtr)/rhsq2;
    double v = (ql[2]/rhsqrtl  + qr[2]/rhsqrtr)/rhsq2;
    double w = (ql[3]/rhsqrtl  + qr[3]/rhsqrtr)/rhsq2;

    double q2 = u*u + v*v + w*w;
// Roe averaged enthalpy
    double enth = ((ql[4]+pl)/rhsqrtl + (qr[4]+pr)/rhsqrtr)/rhsq2;
// speed of sound
    double aa2 = gas_gamma1*(enth - 0.5*q2);
    if (correct && (aa2 < 0))
      aa2 = gas_gamma*minPressure/minDensity;
    double a = std::sqrt(aa2);
// other quantities
    double g1a2 = gas_gamma1/aa2;
    double euv = enth - q2;

// project onto left eigenvectors with Roe averaged values (see Tech Note 1007)
    double a3 = g1a2*(euv*jump[0] + u*jump[1] + v*jump[2] + w*jump[3] - jump[4]);
    double a1 = -v*jump[0] + jump[2];
    double a2 = -w*jump[0] + jump[3];
    double a4 = (jump[1] + (a-u)*jump[0] - a*a3)/(2*a);
    double a0 = jump[0] - a3 - a4;

// compute waves (see Tech Note 1007)

// wave 1: eigenvalue is u-c
    waves(0,0) = a0;
    waves(1,0) = a0*(u-a);
    waves(2,0) = a0*v;
    waves(3,0) = a0*w;
    waves(4,0) = a0*(enth-u*a);
    s[0] = u-a;

// wave 2: 3 eigenvectors for the repeated roots are combined into a
// single wave
    waves(0,1) = a3;
    waves(1,1) = a3*u;
    waves(2,1) = a1 + v*a3;
    waves(3,1) = a2 + w*a3;
    waves(4,1) = a1*v + a2*w + 0.5*a3*q2;
    s[1] = u;

// wave 3: eigenvalue is u+c
    waves(0,2) = a4;
    waves(1,2) = a4*(u+a);
    waves(2,2) = a4*v;
    waves(3,2) = a4*w;
    waves(4,2) = a4*(enth + u*a);
    s[2] = u+a;
  }

  void
  EulerEquation::wavesLax(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    for (unsigned i=0; i<5; ++i)
      waves(i,0) = qr[i]-ql[i];

    double sl[2], sr[2];
    speeds(c, &ql[0], sl);
    speeds(c, &qr[0], sr);
    s[0] = 0.5*(sl[1]+sr[1]); // NOT SURE OF THIS
  }

  double
  EulerEquation::numericalFlux(const Lucee::RectCoordSys& c,
    const double* ql, const double* qr, 
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    double* f)
  {
// NOTE: This numerical flux is using Lax-Fluxes

    double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));

    double fl[5], fr[5];
    flux(c, ql, auxVarsl, fl);
    flux(c, qr, auxVarsr, fr);

    for (unsigned i=0; i<5; ++i)
      f[i] = 0.5*(fr[i]+fl[i]) - 0.5*absMaxs*(qr[i]-ql[i]);

    return absMaxs;
  }

  void
  EulerEquation::qFluctuations(const Lucee::RectCoordSys& c,
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
      double fl[5], fr[5];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));
      double amdqL[5], apdqL[5];
      for (unsigned m=0; m<5; ++m)
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
      
  double
  EulerEquation::pressure(const double* q) const
  {
    double rho1 = 1/getSafeRho(q[0]);
    double pr = (gas_gamma-1)*(q[4] - 0.5*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3])*rho1);
    if (correct && (pr < minPressure))
      return minPressure;
    return pr;
  }

  double
  EulerEquation::getSafeRho(double rho) const
  {
    if (correct && (rho<minDensity))
      return minDensity;
    return rho;
  }
}
