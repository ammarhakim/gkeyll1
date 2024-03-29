/**
 * @file	LcTenMomentEquation.cpp
 *
 * @brief	TenMoment equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTenMomentEquation.h>

// std includes
#include <cmath>

namespace Lucee
{
// set ids for creators
  const char *TenMomentEquation::id = "TenMoment";

// makes indexing easier
  static const unsigned RHO = 0;
  static const unsigned U1 = 1;
  static const unsigned U2 = 2;
  static const unsigned U3 = 3;
  static const unsigned P11 = 4;
  static const unsigned P12 = 5;
  static const unsigned P13 = 6;
  static const unsigned P22 = 7;
  static const unsigned P23 = 8;
  static const unsigned P33 = 9;

  static
  void multiplyByPhiPrime(double p0, double u1, double u2, double u3,
    double p11, double p12, double p13, double p22, double p23, double p33,
    const double wv[10],
    int waveNum, Lucee::Matrix<double>& waves)
  {
    waves(0,waveNum) = wv[0];
    waves(1,waveNum) = p0*wv[1]+u1*wv[0];
    waves(2,waveNum) = p0*wv[2]+u2*wv[0];
    waves(3,waveNum) = p0*wv[3]+u3*wv[0];
    waves(4,waveNum) = wv[4]+2*p0*u1*wv[1]+u1*u1*wv[0];
    waves(5,waveNum) = wv[5]+p0*u1*wv[2]+p0*u2*wv[1]+u1*u2*wv[0];
    waves(6,waveNum) = wv[6]+p0*u1*wv[3]+p0*u3*wv[1]+u1*u3*wv[0];
    waves(7,waveNum) = wv[7]+2*p0*u2*wv[2]+u2*u2*wv[0];
    waves(8,waveNum) = wv[8]+p0*u2*wv[3]+p0*u3*wv[2]+u2*u3*wv[0];
    waves(9,waveNum) = wv[9]+2*p0*u3*wv[3]+u3*u3*wv[0];
  }

  TenMomentEquation::TenMomentEquation()
    : Lucee::HyperEquation(10, 5)
  {
// equation system has 10 equations and 5 waves
  }
  
  void
  TenMomentEquation::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::HyperEquation::readInput(tbl);

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
        Lucee::Except lce("TenMomentEquation::readInput: 'numericalFlux' ");
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
  TenMomentEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToLocal(&inQ[1], &outQ[1]); // momentum
    c.rotateSymMatrixToLocal(&inQ[4], &outQ[4]); // pressure tensor
  }

  void
  TenMomentEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double *inQ, double *outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
    c.rotateSymMatrixToGlobal(&inQ[4], &outQ[4]); // pressure tensor
  }

  void
  TenMomentEquation::flux(const Lucee::RectCoordSys& c, const double* q, 
    const std::vector<const double*>& auxVars, double* f)
  {
    Lucee::FieldPtr<double> v(10);
// compute primitive variables first
    primitive(q, v);
// density flux
    f[RHO] = v[RHO]*v[U1];
// momentum density flux
    f[U1] = v[RHO]*v[U1]*v[U1] + v[P11];
    f[U2] = v[RHO]*v[U1]*v[U2] + v[P12];
    f[U3] = v[RHO]*v[U1]*v[U3] + v[P13];
// total pressure flux
    f[P11] = v[RHO]*v[U1]*v[U1]*v[U1] + 3*v[U1]*v[P11];
    f[P12] = v[RHO]*v[U1]*v[U1]*v[U2] + 2*v[U1]*v[P12] + v[U2]*v[P11];
    f[P13] = v[RHO]*v[U1]*v[U1]*v[U3] + 2*v[U1]*v[P13] + v[U3]*v[P11];
    f[P22] = v[RHO]*v[U1]*v[U2]*v[U2] + v[U1]*v[P22] + 2*v[U2]*v[P12];
    f[P23] = v[RHO]*v[U1]*v[U2]*v[U3] + v[U1]*v[P23] + v[U2]*v[P13] + v[U3]*v[P12];
    f[P33] = v[RHO]*v[U1]*v[U3]*v[U3] + v[U1]*v[P33] + 2*v[U3]*v[P13];
  }

  void
  TenMomentEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s[2])
  {
    double rho = q[RHO];
    double u1 = q[U1]/rho;
    double p11 = q[P11] - rho*u1*u1;
    s[0] = u1-sqrt(3*p11/rho);
    s[1] = u1+sqrt(3*p11/rho);
  }

  double
  TenMomentEquation::maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q)
  {
    double rho = q[RHO];
    double u1 = q[U1]/rho;
    double p11 = q[P11] - rho*u1*u1;
    return std::fabs(u1)+sqrt(3*p11/rho);
  }

  void
  TenMomentEquation::primitive(const double* q, double* v) const
  {
    double rho = q[RHO];
    v[RHO] = rho; // rho
// velocity components
    v[U1] = q[U1]/rho;
    v[U2] = q[U2]/rho;
    v[U3] = q[U3]/rho;
// pressure tensor
    v[P11] = q[P11] - rho*v[U1]*v[U1];
    v[P12] = q[P12] - rho*v[U1]*v[U2];
    v[P13] = q[P13] - rho*v[U1]*v[U3];
    v[P22] = q[P22] - rho*v[U2]*v[U2];
    v[P23] = q[P23] - rho*v[U2]*v[U3];
    v[P33] = q[P33] - rho*v[U3]*v[U3];
  }

  void
  TenMomentEquation::conserved(const double* v, double* q) const
  {
    double rho = v[RHO];
    q[RHO] = rho; // rho
// velocity components
    q[U1] = rho*v[U1];
    q[U2] = rho*v[U2];
    q[U3] = rho*v[U3];
// pressure tensor
    q[P11] = v[P11] + rho*v[U1]*v[U1];
    q[P12] = v[P12] + rho*v[U1]*v[U2];
    q[P13] = v[P13] + rho*v[U1]*v[U3];
    q[P22] = v[P22] + rho*v[U2]*v[U2];
    q[P23] = v[P23] + rho*v[U2]*v[U3];
    q[P33] = v[P33] + rho*v[U3]*v[U3];
  }

  void
  TenMomentEquation::waves(const Lucee::RectCoordSys& c,
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
  TenMomentEquation::wavesRoe(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
// compute right and left primitive values
    double vl[10], vr[10];
    primitive(&ql[0], vl);
    primitive(&qr[0], vr);

// compute averages for use in Riemann solver (Note: These are
// straight averages and not Roe averages. This should not be a
// problem, but might be best to use f-waves to ensure conservation).

    double p0 = 0.5*(vl[RHO]+vr[RHO]);
    double u1 = 0.5*(vl[U1]+vr[U1]);
    double u2 = 0.5*(vl[U2]+vr[U2]);
    double u3 = 0.5*(vl[U3]+vr[U3]);
    double p11 = 0.5*(vl[P11]+vr[P11]);
    double p12 = 0.5*(vl[P12]+vr[P12]);
    double p13 = 0.5*(vl[P13]+vr[P13]);
    double p22 = 0.5*(vl[P22]+vr[P22]);
    double p23 = 0.5*(vl[P23]+vr[P23]);
    double p33 = 0.5*(vl[P33]+vr[P33]);

// The following expressions are cut-paste from the script
// tenmom-eig.mac and then massaged around a bit. Also see Tech-Note
// 1013.

// multiply jumps by (phi')^-1 matrix first (this is because the left
// eigenvectors are computed from the quasilinear form and not the
// primitive form of the equations)
    double phiJump[10];

    phiJump[0] = jump[0];
    phiJump[1] = -(u1*jump[0]-jump[1])/p0;
    phiJump[2] = -(u2*jump[0]-jump[2])/p0;
    phiJump[3] = -(u3*jump[0]-jump[3])/p0;
    phiJump[4] = jump[4]-2*u1*jump[1]+u1*u1*jump[0];
    phiJump[5] = jump[5]-u1*jump[2]+u2*(u1*jump[0]-jump[1]);
    phiJump[6] = jump[6]-u1*jump[3]+u3*(u1*jump[0]-jump[1]);
    phiJump[7] = jump[7]-2*u2*jump[2]+u2*u2*jump[0];
    phiJump[8] = jump[8]-u2*jump[3]+u3*(u2*jump[0]-jump[2]);
    phiJump[9] = jump[9]-2*u3*jump[3]+u3*u3*jump[0];

// bunch of useful stuff
    double sqp0 = sqrt(p0);
    double sqp11 = sqrt(p11);
    double sp11 = p11*p11;
    double fp11 = pow(3,-1.5)*pow(p11,-2.5); // 3**((-3.0)/2.0)*p11**((-5.0)/2.0)
    double sp12 = p12*p12;
    double sp13 = p13*p13;
    double sq3 = sqrt(3.0);
    double thp11 = pow(p11, 1.5);

    double lp[10];
// project jumps on left eigenvectors
    lp[0] = (sqp0*sqp11*(phiJump[1]*p12-phiJump[2]*p11)-phiJump[4]*p12+phiJump[5]*p11)/sp11/2.0;
    lp[1] = (sqp0*sqp11*(phiJump[1]*p13-phiJump[3]*p11)-phiJump[4]*p13+phiJump[6]*p11)/sp11/2.0;
    lp[2] = -(sqp0*sqp11*(phiJump[1]*p12-phiJump[2]*p11)+phiJump[4]*p12-phiJump[5]*p11)/sp11/2.0;
    lp[3] = -(sqp0*sqp11*(phiJump[1]*p13-phiJump[3]*p11)+phiJump[4]*p13-phiJump[6]*p11)/sp11/2.0;
    lp[4] = fp11*(sq3*phiJump[4]*sqp11-3*phiJump[1]*sqp0*p11)/2.0;
    lp[5] = fp11*(3*phiJump[1]*sqp0*p11+sq3*phiJump[4]*sqp11)/2.0;
    lp[6] = (3*phiJump[0]*p11-phiJump[4]*p0)/p11/3.0;
    lp[7] = -(phiJump[4]*p11*p22-4*phiJump[4]*sp12+6*phiJump[5]*p11*p12-3*phiJump[7]*sp11)/sp11/3.0;
    lp[8] = -(phiJump[4]*p11*p23+(3*phiJump[5]*p11-4*phiJump[4]*p12)*p13+3*phiJump[6]*p11*p12-3*phiJump[8]*sp11)/sp11/3.0;
    lp[9] = -(phiJump[4]*p11*p33-4*phiJump[4]*sp13+6*phiJump[6]*p11*p13-3*phiJump[9]*sp11)/sp11/3.0;

// compute waves
    double wv[10];

// wave 1: eigenvalues 1,2
    s[0] = u1-sqrt(p11/p0);

    wv[0] = 0;
    wv[1] = 0;
    wv[2] = -lp[0]*sqp11/sqp0;
    wv[3] = -lp[1]*sqp11/sqp0;
    wv[4] = 0;
    wv[5] = lp[0]*p11;
    wv[6] = lp[1]*p11;
    wv[7] = 2*lp[0]*p12;
    wv[8] = lp[0]*p13+lp[1]*p12;
    wv[9] = 2*lp[1]*p13;

    multiplyByPhiPrime(p0, u1, u2, u3, p11, p12, p13, p22, p23, p33, wv, 0, waves);

// wave 2: eigenvalues 3,4
    s[1] = u1+sqrt(p11/p0);

    wv[0] = 0;
    wv[1] = 0;
    wv[2] = lp[2]*sqp11/sqp0;
    wv[3] = lp[3]*sqp11/sqp0;
    wv[4] = 0;
    wv[5] = lp[2]*p11;
    wv[6] = lp[3]*p11;
    wv[7] = 2*lp[2]*p12;
    wv[8] = lp[2]*p13+lp[3]*p12;
    wv[9] = 2*lp[3]*p13;

    multiplyByPhiPrime(p0, u1, u2, u3, p11, p12, p13, p22, p23, p33, wv, 1, waves);

// wave 3: eigenvalue 5
    s[2] = u1-sqrt(3*p11/p0);

    wv[0] = lp[4]*p0*p11;
    wv[1] = -sq3*lp[4]*thp11/sqp0;
    wv[2] = -sq3*lp[4]*sqp11*p12/sqp0;
    wv[3] = -sq3*lp[4]*sqp11*p13/sqp0;
    wv[4] = 3*lp[4]*sp11;
    wv[5] = 3*lp[4]*p11*p12;
    wv[6] = 3*lp[4]*p11*p13;
    wv[7] = lp[4]*(p11*p22+2*sp12);
    wv[8] = lp[4]*(p11*p23+2*p12*p13);
    wv[9] = lp[4]*(p11*p33+2*sp13);

    multiplyByPhiPrime(p0, u1, u2, u3, p11, p12, p13, p22, p23, p33, wv, 2, waves);

// wave 4: eigenvalue 6
    s[3] = u1+sqrt(3*p11/p0);

    wv[0] = lp[5]*p0*p11;
    wv[1] = sq3*lp[5]*thp11/sqrt(p0);
    wv[2] = sq3*lp[5]*sqp11*p12/sqp0;
    wv[3] = sq3*lp[5]*sqp11*p13/sqp0;
    wv[4] = 3*lp[5]*sp11;
    wv[5] = 3*lp[5]*p11*p12;
    wv[6] = 3*lp[5]*p11*p13;
    wv[7] = lp[5]*(p11*p22+2*sp12);
    wv[8] = lp[5]*(p11*p23+2*p12*p13);
    wv[9] = lp[5]*(p11*p33+2*sp13);

    multiplyByPhiPrime(p0, u1, u2, u3, p11, p12, p13, p22, p23, p33, wv, 3, waves);

// wave 5: eigenvalue 7,8,9,10
    s[4] = u1;

    wv[0] = lp[6];
    wv[1] = 0;
    wv[2] = 0;
    wv[3] = 0;
    wv[4] = 0;
    wv[5] = 0;
    wv[6] = 0;
    wv[7] = lp[7];
    wv[8] = lp[8];
    wv[9] = lp[9];

    multiplyByPhiPrime(p0, u1, u2, u3, p11, p12, p13, p22, p23, p33, wv, 4, waves);
  }

  void
  TenMomentEquation::wavesLax(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    if (useIntermediateWave)
    {
// this uses the HLLE technique to construct an intermediate state
      std::vector<const double*> auxVars;
      double fl[10], fr[10];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[0]+sr[0]);
      s[1] = 0.5*(sl[1]+sr[1]);

// compute intermediate HLLE state
      double sdiff1 = 1/(s[1]-s[0]);
      double qHHLE[10];
      for (unsigned i=0; i<10; ++i)
        qHHLE[i] = (s[1]*qr[i]-s[0]*ql[i]+fl[i]-fr[i])*sdiff1;

// compute waves
      for (unsigned i=0; i<10; ++i)
      {
        waves(i,0) = qHHLE[i]-ql[i];
        waves(i,1) = qr[i]-qHHLE[i];
      }
    }
    else
    {
// use a single wave
      for (unsigned i=0; i<10; ++i)
        waves(i,0) = qr[i]-ql[i];
      
      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[1]+sr[1]);
    }
  }

  void
  TenMomentEquation::qFluctuations(const Lucee::RectCoordSys& c,
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
      double fl[10], fr[10];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));
      double amdqL[10], apdqL[10];
      for (unsigned m=0; m<10; ++m)
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
  TenMomentEquation::isInvariantDomain(const double* q) const
  {
    double rho = q[RHO];
    if (rho<=0.0)
      return false;

    double u = q[U1]/rho;
    double v = q[U2]/rho;
    double w = q[U3]/rho;

    if ((q[P11] - rho*u*u)<=0)
      return false;
    if ((q[P22] - rho*v*v)<=0)
      return false;
    if ((q[P33] - rho*w*w)<=0)
      return false;

    return true;
  }
}
