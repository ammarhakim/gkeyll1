/**
 * @file	LcTwentyMomentEquation.cpp
 *
 * @brief	TwentyMoment equations for gas-dynamics.
 */
// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcTwentyMomentEquation.h>
#include <LcLinAlgebra.h>
#include <LcMatrix.h>

// std includes
#include <cmath>

// eigen includes
#include <Eigen/Eigen>

namespace Lucee
{

// set ids for creators
  const char *TwentyMomentEquation::id = "TwentyMoment";
    
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
  static const unsigned Q111 = 10;
  static const unsigned Q112 = 11;
  static const unsigned Q113 = 12;
  static const unsigned Q122 = 13;
  static const unsigned Q123 = 14;
  static const unsigned Q133 = 15;
  static const unsigned Q222 = 16;
  static const unsigned Q223 = 17;
  static const unsigned Q233 = 18;
  static const unsigned Q333 = 19;
  
  static
  void multiplyByPhiPrime(double p0, double u1, double u2, double u3,
    double p11, double p12, double p13, double p22, double p23, double p33,
    const double wv[20],
    int waveNum, Lucee::Matrix<double>& waves)
  {
    waves(0,waveNum) = wv[0];
    waves(1,waveNum) = u1*wv[0] + p0*wv[1];
    waves(2,waveNum) = u2*wv[0] + p0*wv[2];
    waves(3,waveNum) = u3*wv[0] + p0*wv[3];
    waves(4,waveNum) = u1*u1*wv[0] + 2*p0*u1*wv[1] + wv[4];
    waves(5,waveNum) = u1*u2*wv[0] + p0*u2*wv[1] + p0*u1*wv[2] + wv[5];
    waves(6,waveNum) = u1*u3*wv[0] + p0*u3*wv[1] + p0*u1*wv[3] + wv[6];
    waves(7,waveNum) = u2*u2*wv[0] + 2*p0*u2*wv[2] + wv[7];
    waves(8,waveNum) = u2*u3*wv[0] + p0*u3*wv[2] +  p0*u2*wv[3] + wv[8];
    waves(9,waveNum) = u3*u3*wv[0] + 2*p0*u3*wv[3] + wv[9];
    waves(10,waveNum) = u1*u1*u1*wv[0] + (3*p11 + 3*p0*u1*u1)*wv[1] + 3*u1*wv[4] + wv[10];
    waves(11,waveNum) = u1*u1*u2*wv[0] + (2*p12 + 2*p0*u1*u2)*wv[1] +  (p11 + p0*u1*u1)*wv[2] + u2*wv[4] + 2*u1*wv[5] + wv[11];
    waves(12,waveNum) = u1*u1*u3*wv[0] + (2*p13 + 2*p0*u1*u3)*wv[1] +  (p11 + p0*u1*u1)*wv[3] + u3*wv[4] + 2*u1*wv[6] + wv[12];
    waves(13,waveNum) = u1*u2*u2*wv[0] + (p22 + p0*u2*u2)*wv[1] +  (2*p12 + 2*p0*u1*u2)*wv[2] + 2*u2*wv[5] + u1*wv[7] + wv[13];
    waves(14,waveNum) = u1*u2*u3*wv[0] + (p23 + p0*u2*u3)*wv[1] +  (p13 + p0*u1*u3)*wv[2] + (p12 + p0*u1*u2)*wv[3] + u3*wv[5] +  u2*wv[6] + u1*wv[8] + wv[14];
    waves(15,waveNum) = u1*u3*u3*wv[0] +  (p33 + p0*u3*u3)*wv[1] + (2*p13 + 2*p0*u1*u3)*wv[3] +  2*u3*wv[6] + u1*wv[9] + wv[15];
    waves(16,waveNum) = u2*u2*u2*wv[0] + (3*p22 + 3*p0*u2*u2)*wv[2] + 3*u2*wv[7] + wv[16];
    waves(17,waveNum) = u2*u2*u3*wv[0] + (2*p23 + 2*p0*u2*u3)*wv[2] +  (p22 + p0*u2*u2)*wv[3] + u3*wv[7] + 2*u2*wv[8] + wv[17];
    waves(18,waveNum) = u2*u3*u3*wv[0] + (p33 + p0*u3*u3)*wv[2] +  (2*p23 + 2*p0*u2*u3)*wv[3] + 2*u3*wv[8] + u2*wv[9] + wv[18];
    waves(19,waveNum) = u3*u3*u3*wv[0] + (3*p33 + 3*p0*u3*u3)*wv[3] + 3*u3*wv[9] + wv[19];
  }

  TwentyMomentEquation::TwentyMomentEquation()
    : Lucee::HyperEquation(20,20)
    {
// twenty equations, 9 waves
    }

  void TwentyMomentEquation::readInput(Lucee::LuaTable& tbl)
  {
    // call base class method
    Lucee::HyperEquation::readInput(tbl);

    numFlux = NF_ROE;
    if(tbl.hasString("numericalFlux")) {
        std::string nf = tbl.getString("numericalFlux");
        if(nf == "roe")
            numFlux = NF_ROE;
        else if(nf == "lax")
            numFlux = NF_LAX;
        else {
            Lucee::Except lce("TwentyMomentEquation::readInput: 'numericalFlux' ");
            lce << nf << " not recognized!" << std::endl;
            throw lce;
        }
    }

    if(numFlux == NF_LAX) {
        useIntermediateWave = false;
        if(tbl.hasBool("useIntermediateWave"))
            // this flag activates the intermediate state wave: it is best to use
            // this when using Lax fluxes with second order wave-propagation
            // scheme so as to avoid asymmetries.
            useIntermediateWave = tbl.getBool("useIntermediateWave");

        // adjust number of waves accordingly
        if(useIntermediateWave)
            this->setNumWaves(2);
        else
            this->setNumWaves(1);
    }
    
  }

  void TwentyMomentEquation::rotateToLocal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ)
  {
    outQ[0] = inQ[0];
    c.rotateVecToLocal(&inQ[1], &outQ[1]);       // momentum
    c.rotateSymMatrixToLocal(&inQ[4], &outQ[4]); // pressure tensor
    c.rotateSymTensorToLocal(&inQ[10], &outQ[10]); // heat flux tensor
  }


  void TwentyMomentEquation::rotateToGlobal(const Lucee::RectCoordSys& c, const double* inQ, double* outQ)
  {
    outQ[0] = inQ[0]; // density
    c.rotateVecToGlobal(&inQ[1], &outQ[1]); // momentum
    c.rotateSymMatrixToGlobal(&inQ[4], &outQ[4]); // pressure tensor
    c.rotateSymTensorToGlobal(&inQ[10], &outQ[10]); // heat flux tensor
  }


  void TwentyMomentEquation::flux(const Lucee::RectCoordSys& c, const double* q, const std::vector<const double*>& auxVars, double* f)
  {
      Lucee::FieldPtr<double> v(20);
      // compute primitive variables first
      primitive(q, v);
      // density flux
      f[RHO] = v[RHO] * v[U1];
      // momentum density flux
      f[U1] = v[RHO] * v[U1] * v[U1] + v[P11];
      f[U2] = v[RHO] * v[U1] * v[U2] + v[P12];
      f[U3] = v[RHO] * v[U1] * v[U3] + v[P13];
      // total pressure flux
      f[P11] = v[RHO] * v[U1] * v[U1] * v[U1] + 3 * v[U1] * v[P11] + v[Q111];
      f[P12] = v[RHO] * v[U1] * v[U1] * v[U2] + 2 * v[U1] * v[P12] + v[U2] * v[P11] + v[Q112];
      f[P13] = v[RHO] * v[U1] * v[U1] * v[U3] + 2 * v[U1] * v[P13] + v[U3] * v[P11] + v[Q113];
      f[P22] = v[RHO] * v[U1] * v[U2] * v[U2] + v[U1] * v[P22] + 2 * v[U2] * v[P12] + v[Q122];
      f[P23] = v[RHO] * v[U1] * v[U2] * v[U3] + v[U1] * v[P23] + v[U2] * v[P13] + v[U3] * v[P12] + v[Q123];
      f[P33] = v[RHO] * v[U1] * v[U3] * v[U3] + v[U1] * v[P33] + 2 * v[U3] * v[P13] + v[Q133];
      // make the inverse of the pressure tensor
      /*      double det = -v[P13]*v[P13]*v[P22] + 2*v[P12]*v[P13]*v[P23] -v[P11]*v[P23]*v[P23] - v[P12]*v[P12]*v[P33] + v[P11]*v[P22]*v[P33];
      double PInv[6];
      PInv[0] = (-v[P23]*v[P23] + v[P22]*v[P33])/det;
      PInv[1] = (v[P13]*v[P23] - v[P12]*v[P33])/det;
      PInv[2] = (-v[P13]*v[P22] + v[P12]*v[P23])/det;
      PInv[3] = (-v[P13]*v[P13] + v[P11]*v[P33])/det;
      PInv[4] = (v[P12]*v[P13] - v[P11]*v[P23])/det;
      PInv[5] = (-v[P12]*v[P12] + v[P11]*v[P22])/det; */
      // total heat flux flux
      double k; // parameter which measures deviation from maxwellian (sort of Levermore, sort of Pearson). Totally not rigourous 
      k = (1 + (v[RHO]*(v[Q111]*v[Q111]/(v[P11]*v[P11]*v[P11]) + v[Q222]*v[Q222]/(v[P22]*v[P22]*v[P22]) + v[Q333]*v[Q333]/(v[P33]*v[P33]*v[P33])))/2 );
      f[Q111] = 4*v[Q111]*v[U1] + (3*k*(v[P11]*v[P11]))/(v[RHO]) + 6*v[P11]*(v[U1]*v[U1]) + v[RHO]*(v[U1]*v[U1]*v[U1]*v[U1]);

      f[Q112] = (3*k*v[P11]*v[P12])/(v[RHO]) + 3*v[Q112]*v[U1] + v[Q111]*v[U2] + 3*v[P11]*v[U1]*v[U2] + 3*v[P12]*(v[U1]*v[U1]) + v[RHO]*v[U2]*(v[U1]*v[U1]*v[U1]);

      f[Q113] = (3*k*v[P11]*v[P13])/(v[RHO]) + 3*v[Q113]*v[U1] + v[Q111]*v[U3] + 3*v[P11]*v[U1]*v[U3] + 3*v[P13]*(v[U1]*v[U1]) + v[RHO]*v[U3]*(v[U1]*v[U1]*v[U1]);

      f[Q122] = 2*v[Q122]*v[U1] + 2*v[Q112]*v[U2] + 4*v[P12]*v[U1]*v[U2] + k*(v[P11]*v[P22] + 2*(v[P12]*v[P12]))/(v[RHO]) + v[P22]*(v[U1]*v[U1]) + v[P11]*(v[U2]*v[U2]) + v[RHO]*(v[U1]*v[U1])*(v[U2]*v[U2]);

      f[Q123] = (2*v[P12]*(k*v[P13] + v[RHO]*v[U1]*v[U3]) + v[P11]*(k*v[P23] + v[RHO]*v[U2]*v[U3]) + v[RHO]*(2*v[Q123]*v[U1] + v[Q113]*v[U2] + 2*v[P13]*v[U1]*v[U2] + v[Q112]*v[U3] + v[P23]*(v[U1]*v[U1]) + v[RHO]*v[U2]*v[U3]*(v[U1]*v[U1])))/(v[RHO]);

      f[Q133] = 2*v[Q133]*v[U1] + 2*v[Q113]*v[U3] + 4*v[P13]*v[U1]*v[U3] + k*(v[P11]*v[P33] + 2*(v[P13]*v[P13]))/(v[RHO]) + v[P33]*(v[U1]*v[U1]) + v[P11]*(v[U3]*v[U3]) + v[RHO]*(v[U1]*v[U1])*(v[U3]*v[U3]);

      f[Q222] = (3*k*v[P12]*v[P22])/(v[RHO]) + v[Q222]*v[U1] + 3*v[U2]*(v[Q122] + v[P22]*v[U1] + v[P12]*v[U2]) + v[RHO]*v[U1]*(v[U2]*v[U2]*v[U2]);

      f[Q223] = (2*v[P12]*(k*v[P23] + v[RHO]*v[U2]*v[U3]) + v[P13]*(k*v[P22] + v[RHO]*(v[U2]*v[U2])) + v[RHO]*(v[Q223]*v[U1] + 2*v[Q123]*v[U2] + 2*v[P23]*v[U1]*v[U2] + v[Q122]*v[U3] + v[P22]*v[U1]*v[U3] + v[RHO]*v[U1]*v[U3]*(v[U2]*v[U2])))/(v[RHO]);

      f[Q233] = (2*v[P13]*(k*v[P23] + v[RHO]*v[U2]*v[U3]) + v[P12]*(k*v[P33] + v[RHO]*(v[U3]*v[U3])) + v[RHO]*(v[Q233]*v[U1] + v[Q133]*v[U2] + v[P33]*v[U1]*v[U2] + 2*v[Q123]*v[U3] + 2*v[P23]*v[U1]*v[U3] + v[RHO]*v[U1]*v[U2]*(v[U3]*v[U3])))/(v[RHO]);

      f[Q333] = (3*k*v[P13]*v[P33])/(v[RHO]) + v[Q333]*v[U1] + 3*v[U3]*(v[Q133] + v[P33]*v[U1] + v[P13]*v[U3]) + v[RHO]*v[U1]*(v[U3]*v[U3]*v[U3]);
  }

  void TwentyMomentEquation::speeds(const Lucee::RectCoordSys& c, const double* q, double s [2])
  {
    double rho = q[RHO];
    double u1 = q[U1]/rho;
    double p11 = q[P11] - rho*u1*u1;
    double v = sqrt((3 + sqrt(6))*p11/rho)*1.5;
    s[0] = u1 - v;
    s[1] = u1 + v;
  }

  double TwentyMomentEquation::maxAbsSpeed(const Lucee::RectCoordSys& c, const double* q)
  {    
    double rho = q[RHO];
    double u1 = q[U1]/rho;
    double p11 = q[P11] - rho*u1*u1;
    return std::fabs(u1) + sqrt((3 + sqrt(6))*p11/rho)*1.5; 
  }

  void TwentyMomentEquation::primitive(const double* q, double* v) const
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
// heat flux tensor
    v[Q111] = q[Q111] - 3*v[P11]*v[U1] - rho*v[U1]*v[U1]*v[U1];
    v[Q112] = q[Q112] - 2*v[P12]*v[U1] - v[P11]*v[U2] - rho*v[U1]*v[U1]*v[U2];
    v[Q113] = q[Q113] - 2*v[P13]*v[U1] - v[P11]*v[U3] - rho*v[U1]*v[U1]*v[U3];
    v[Q122] = q[Q122] - v[P22]*v[U1] - 2*v[P12]*v[U2] - rho*v[U1]*v[U2]*v[U2];
    v[Q123] = q[Q123] - v[P23]*v[U1] - v[P13]*v[U2] - v[P12]*v[U3] - rho*v[U1]*v[U2]*v[U3];
    v[Q133] = q[Q133] - v[P33]*v[U1] - 2*v[P13]*v[U3] - rho*v[U1]*v[U3]*v[U3];
    v[Q222] = q[Q222] - 3*v[P22]*v[U2] - rho*v[U2]*v[U2]*v[U2];
    v[Q223] = q[Q223] - 2*v[P23]*v[U2] - v[P22]*v[U3] - rho*v[U2]*v[U2]*v[U3];
    v[Q233] = q[Q233] - v[P33]*v[U2] - 2*v[P23]*v[U3] - rho*v[U2]*v[U3]*v[U3];
    v[Q333] = q[Q333] - 3*v[P33]*v[U3] - rho*v[U3]*v[U3]*v[U3];
  }

  void TwentyMomentEquation::conserved(const double* v, double* q) const
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
// heat flux tensor
    q[Q111] = v[Q111] + 3*v[P11]*v[U1] + rho*v[U1]*v[U1]*v[U1];
    q[Q112] = v[Q112] + 2*v[P12]*v[U1] + v[P11]*v[U2] + rho*v[U1]*v[U1]*v[U2];
    q[Q113] = v[Q113] + 2*v[P13]*v[U1] + v[P11]*v[U3] + rho*v[U1]*v[U1]*v[U3];
    q[Q122] = v[Q122] + v[P22]*v[U1] + 2*v[P12]*v[U2] + rho*v[U1]*v[U2]*v[U2];
    q[Q123] = v[Q123] + v[P23]*v[U1] + v[P13]*v[U2] + v[P12]*v[U3] + rho*v[U1]*v[U2]*v[U3];
    q[Q133] = v[Q133] + v[P33]*v[U1] + 2*v[P13]*v[U3] + rho*v[U1]*v[U3]*v[U3];
    q[Q222] = v[Q222] + 3*v[P22]*v[U2] + rho*v[U2]*v[U2]*v[U2];
    q[Q223] = v[Q223] + 2*v[P23]*v[U2] + v[P22]*v[U3] + rho*v[U2]*v[U2]*v[U3];
    q[Q233] = v[Q233] + v[P33]*v[U2] + 2*v[P23]*v[U3] + rho*v[U2]*v[U3]*v[U3];
    q[Q333] = v[Q333] + 3*v[P33]*v[U3] + rho*v[U3]*v[U3]*v[U3];
  }


  void TwentyMomentEquation::waves(const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& jump,
    const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr,
    const std::vector<const double*>& auxVarsl, const std::vector<const double*>& auxVarsr,
    Lucee::Matrix<double>& waves,    
    Lucee::FieldPtr<double>& s)
  {
    if (numFlux == NF_ROE)
      wavesRoe(c, jump, ql, qr, waves, s);
    else if (numFlux == NF_LAX)
      wavesLax(c, jump, ql, qr, waves, s);
    else
    { /* this can't happen */ }
  }

  void TwentyMomentEquation::wavesRoe(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& jump, const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr, Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
// compute right and left primitive values
    double vl[20], vr[20];
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
    double q111 = 0.5*(vl[Q111]+vr[Q111]);
    double q112 = 0.5*(vl[Q112]+vr[Q112]);
    double q113 = 0.5*(vl[Q113]+vr[Q113]);
    double q122 = 0.5*(vl[Q122]+vr[Q122]);
    double q123 = 0.5*(vl[Q123]+vr[Q123]);
    double q133 = 0.5*(vl[Q133]+vr[Q133]);
    double q222 = 0.5*(vl[Q222]+vr[Q222]);
    double q223 = 0.5*(vl[Q223]+vr[Q223]);
    double q233 = 0.5*(vl[Q233]+vr[Q233]);
    double q333 = 0.5*(vl[Q333]+vr[Q333]);


// multiply jumps by (phi')^-1 matrix first (this is because the left
// eigenvectors are computed from the quasilinear form and not the
// primitive form of the equations)
    
// this was generated by mathematica and cleaned up.
/*  eigen code which does not seem to work numerically though it solves the same equation

    Eigen::VectorXd phiJump(20);
    phiJump(0) = jump[0];
    phiJump(1) = (-(u1*jump[0]) + jump[1])/p0;
    phiJump(2) = (-(u2*jump[0]) + jump[2])/p0;
    phiJump(3) = (-(u3*jump[0]) + jump[3])/p0;
    phiJump(4) = u1*u1*jump[0] - 2*u1*jump[1] + jump[4];
    phiJump(5) = u1*u2*jump[0] - u2*jump[1] - u1*jump[2] + jump[5];
    phiJump(6) = u1*u3*jump[0] - u3*jump[1] - u1*jump[3] + jump[6];
    phiJump(7) = u2*u2*jump[0] - 2*u2*jump[2] + jump[7];
    phiJump(8) = u2*u3*jump[0] - u3*jump[2] - u2*jump[3] + jump[8];
    phiJump(9) = u3*u3*jump[0] - 2*u3*jump[3] + jump[9];
    phiJump(10) =  -(u1*u1*u1*jump[0]) + (3*p11*(u1*jump[0] - jump[1]))/p0 + 3*u1*u1*jump[1] - 3*u1*jump[4] + jump[10];
    phiJump(11) = (2*p12*(u1*jump[0] - jump[1]) + p11*(u2*jump[0] - jump[2]) + p0*(u1*u1*(-(u2*jump[0]) + jump[2]) - u2*jump[4] + 2*u1*(u2*jump[1] - jump[5]) + jump[11]))/p0;
    phiJump(12) = (2*p13*(u1*jump[0] - jump[1]) + p11*(u3*jump[0] - jump[3]) +    p0*(u1*u1*(-(u3*jump[0]) + jump[3]) - u3*jump[4] + 2*u1*(u3*jump[1] - jump[6]) + jump[12]))/p0;
    phiJump(13) =  (p22*(u1*jump[0] - jump[1]) + 2*p12*(u2*jump[0] - jump[2]) +    p0*(u2*u2*jump[1] - 2*u2*jump[5] - u1*(u2*u2*jump[0] - 2*u2*jump[2] +        jump[7]) + jump[13]))/p0;
    phiJump(14) =  (p23*u1*jump[0] + p13*u2*jump[0] + p12*u3*jump[0] - p0*u1*u2*u3*jump[0] - p23*jump[1] + p0*u2*u3*jump[1] - p13*jump[2] + p0*u1*u3*jump[2] - p12*jump[3] + p0*u1*u2*jump[3] - p0*u3*jump[5] - p0*u2*jump[6] - p0*u1*jump[8] + p0*jump[14])/p0;
    phiJump(15) =  (p33*(u1*jump[0] - jump[1]) + 2*p13*(u3*jump[0] - jump[3]) +    p0*(u3*u3*jump[1] - 2*u3*jump[6] - u1*(u3*u3*jump[0] - 2*u3*jump[3] + jump[9]) + jump[15]))/p0;
    phiJump(16) = -(u2*u2*u2*jump[0]) +   (3*p22*(u2*jump[0] - jump[2]))/p0 + 3*u2*u2*jump[2] - 3*u2*jump[7] + jump[16];
    phiJump(17) = (2*p23*(u2*jump[0] - jump[2]) + p22*(u3*jump[0] - jump[3]) +    p0*(u2*u2*(-(u3*jump[0]) + jump[3]) - u3*jump[7] + 2*u2*(u3*jump[2] - jump[8]) + jump[17]))/p0;
    phiJump(18) =  (p33*(u2*jump[0] - jump[2]) + 2*p23*(u3*jump[0] - jump[3]) +    p0*(u3*u3*jump[2] - 2*u3*jump[8] - u2*(u3*u3*jump[0] - 2*u3*jump[3] + jump[9]) + jump[18]))/p0;
    phiJump(19) = -(u3*u3*u3*jump[0]) + (3*p33*(u3*jump[0] - jump[3]))/p0 + 3*u3*u3*jump[3] - 3*u3*jump[9] + jump[19];
*/

    Lucee::Matrix<double> phiJump(20,1);
    phiJump(0,0) = jump[0];
    phiJump(1,0) = (-(u1*jump[0]) + jump[1])/p0;
    phiJump(2,0) = (-(u2*jump[0]) + jump[2])/p0;
    phiJump(3,0) = (-(u3*jump[0]) + jump[3])/p0;
    phiJump(4,0) = u1*u1*jump[0] - 2*u1*jump[1] + jump[4];
    phiJump(5,0) = u1*u2*jump[0] - u2*jump[1] - u1*jump[2] + jump[5];
    phiJump(6,0) = u1*u3*jump[0] - u3*jump[1] - u1*jump[3] + jump[6];
    phiJump(7,0) = u2*u2*jump[0] - 2*u2*jump[2] + jump[7];
    phiJump(8,0) = u2*u3*jump[0] - u3*jump[2] - u2*jump[3] + jump[8];
    phiJump(9,0) = u3*u3*jump[0] - 2*u3*jump[3] + jump[9];
    phiJump(10,0) =  -(u1*u1*u1*jump[0]) + (3*p11*(u1*jump[0] - jump[1]))/p0 + 3*u1*u1*jump[1] - 3*u1*jump[4] + jump[10];
    phiJump(11,0) = (2*p12*(u1*jump[0] - jump[1]) + p11*(u2*jump[0] - jump[2]) + p0*(u1*u1*(-(u2*jump[0]) + jump[2]) - u2*jump[4] + 2*u1*(u2*jump[1] - jump[5]) + jump[11]))/p0;
    phiJump(12,0) = (2*p13*(u1*jump[0] - jump[1]) + p11*(u3*jump[0] - jump[3]) +    p0*(u1*u1*(-(u3*jump[0]) + jump[3]) - u3*jump[4] + 2*u1*(u3*jump[1] - jump[6]) + jump[12]))/p0;
    phiJump(13,0) =  (p22*(u1*jump[0] - jump[1]) + 2*p12*(u2*jump[0] - jump[2]) +    p0*(u2*u2*jump[1] - 2*u2*jump[5] - u1*(u2*u2*jump[0] - 2*u2*jump[2] +        jump[7]) + jump[13]))/p0;
    phiJump(14,0) =  (p23*u1*jump[0] + p13*u2*jump[0] + p12*u3*jump[0] - p0*u1*u2*u3*jump[0] - p23*jump[1] + p0*u2*u3*jump[1] - p13*jump[2] + p0*u1*u3*jump[2] - p12*jump[3] + p0*u1*u2*jump[3] - p0*u3*jump[5] - p0*u2*jump[6] - p0*u1*jump[8] + p0*jump[14])/p0;
    phiJump(15,0) =  (p33*(u1*jump[0] - jump[1]) + 2*p13*(u3*jump[0] - jump[3]) +    p0*(u3*u3*jump[1] - 2*u3*jump[6] - u1*(u3*u3*jump[0] - 2*u3*jump[3] + jump[9]) + jump[15]))/p0;
    phiJump(16,0) = -(u2*u2*u2*jump[0]) +   (3*p22*(u2*jump[0] - jump[2]))/p0 + 3*u2*u2*jump[2] - 3*u2*jump[7] + jump[16];
    phiJump(17,0) = (2*p23*(u2*jump[0] - jump[2]) + p22*(u3*jump[0] - jump[3]) +    p0*(u2*u2*(-(u3*jump[0]) + jump[3]) - u3*jump[7] + 2*u2*(u3*jump[2] - jump[8]) + jump[17]))/p0;
    phiJump(18,0) =  (p33*(u2*jump[0] - jump[2]) + 2*p23*(u3*jump[0] - jump[3]) +    p0*(u3*u3*jump[2] - 2*u3*jump[8] - u2*(u3*u3*jump[0] - 2*u3*jump[3] + jump[9]) + jump[18]))/p0;
    phiJump(19,0) = -(u3*u3*u3*jump[0]) + (3*p33*(u3*jump[0] - jump[3]))/p0 + 3*u3*u3*jump[3] - 3*u3*jump[9] + jump[19];


    Lucee::Matrix<double> A(20,20); // quasilinear flux matrix
    Lucee::Matrix<double> phiprime(20,20);
    /* eigen code
    Eigen::MatrixXd A; // quasilinear flux matrix
    Eigen::MatrixXd phiprime;
    A = Eigen::MatrixXd::Constant(20,20,0.0);
    phiprime = Eigen::MatrixXd::Constant(20,20,0.0);
    */

    double k; // parameter which measures deviation from maxwellian (sort of Levermore, sort of Pearson). Totally not rigourous 
    // here I will attempt this using k as the heat-flux deviation from a maxwellian.

    k = (1 + (p0*q111*q111/(p11*p11*p11) + p0*q222*q222/(p22*p22*p22) + p0*q333*q333/(p33*p33*p33))/2 );
    /*    if (k > 5) {
      std::cout << "k = " << k << "\n";
      }*/
    A(0,0) = u1; A(0,1) = p0; 
    A(1,4) = 1/p0; A(1,1) = u1; 
    A(2,5) = 1/p0; A(2,2) = u1; 
    A(3,6) = 1/p0; A(3,3) = u1; 
    A(4,4) = u1; A(4,10) = 1; A(4,1) = 3*p11; 
    A(5,5) = u1; A(5,11) = 1; A(5,1) = 2*p12; A(5,2) = p11; 
    A(6,6) = u1; A(6,12) = 1; A(6,1) = 2*p13; A(6,3) = p11; 
    A(7,7) = u1; A(7,13) = 1; A(7,1) = p22; A(7,2) = 2*p12; 
    A(8,8) = u1; A(8,14) = 1; A(8,1) = p23; A(8,2) = p13; A(8,3) = p12; 
    A(9,9) = u1; A(9,15) = 1; A(9,1) = p33; A(9,3) = 2*p13; 
    A(10,0) = (-3*k*(p11*p11))/(p0*p0); A(10,4) = (-3*p11)/p0 + (6*k*p11)/p0; A(10,10) = u1; A(10,1) = 4*q111; 
    A(11,0) = (-3*k*p11*p12)/(p0*p0); A(11,4) = (-2*p12)/p0 + (3*k*p12)/p0; A(11,5) = -(p11/p0) + (3*k*p11)/p0; A(11,11) = u1; A(11,1) = 3*q112; A(11,2) = q111; 
    A(12,0) = (-3*k*p11*p13)/(p0*p0); A(12,4) = (-2*p13)/p0 + (3*k*p13)/p0; A(12,6) = -(p11/p0) + (3*k*p11)/p0; A(12,12) = u1; A(12,1) = 3*q113; A(12,3) = q111; 
    A(13,0) = (-2*k*(p12*p12))/(p0*p0) - (k*p11*p22)/(p0*p0); A(13,4) = -(p22/p0) + (k*p22)/p0; A(13,5) = (-2*p12)/p0 + (4*k*p12)/p0; A(13,7) = (k*p11)/p0; A(13,13) = u1; A(13,1) = 2*q122; A(13,2) = 2*q112; 
    A(14,0) = (-2*k*p12*p13)/(p0*p0) - (k*p11*p23)/(p0*p0); A(14,4) = -(p23/p0) + (k*p23)/p0; A(14,5) = -(p13/p0) + (2*k*p13)/p0; A(14,6) = -(p12/p0) + (2*k*p12)/p0; A(14,8) = (k*p11)/p0; A(14,14) = u1; A(14,1) = 2*q123; A(14,2) = q113; A(14,3) = q112; 
    A(15,0) = (-2*k*(p13*p13))/(p0*p0) - (k*p11*p33)/(p0*p0); A(15,4) = -(p33/p0) + (k*p33)/p0; A(15,6) = (-2*p13)/p0 + (4*k*p13)/p0; A(15,9) = (k*p11)/p0; A(15,15) = u1; A(15,1) = 2*q133; A(15,3) = 2*q113; 
    A(16,0) = (-3*k*p12*p22)/(p0*p0); A(16,5) = (-3*p22)/p0 + (3*k*p22)/p0; A(16,7) = (3*k*p12)/p0; A(16,16) = u1; A(16,1) = q222; A(16,2) = 3*q122; 
    A(17,0) = -((k*p13*p22)/(p0*p0)) - (2*k*p12*p23)/(p0*p0); A(17,5) = (-2*p23)/p0 + (2*k*p23)/p0; A(17,6) = -(p22/p0) + (k*p22)/p0; A(17,7) = (k*p13)/p0; A(17,8) = (2*k*p12)/p0; A(17,17) = u1; A(17,1) = q223; A(17,2) = 2*q123; A(17,3) = q122; 
    A(18,0) = (-2*k*p13*p23)/(p0*p0) - (k*p12*p33)/(p0*p0); A(18,5) = -(p33/p0) + (k*p33)/p0; A(18,6) = (-2*p23)/p0 + (2*k*p23)/p0; A(18,8) = (2*k*p13)/p0; A(18,9) = (k*p12)/p0; A(18,18) = u1; A(18,1) = q233; A(18,2) = q133; A(18,3) = 2*q123; 
    A(19,0) = (-3*k*p13*p33)/(p0*p0); A(19,6) = (-3*p33)/p0 + (3*k*p33)/p0; A(19,9) = (3*k*p13)/p0; A(19,19) = u1; A(19,1) = q333; A(19,3) = 3*q133; 
 

    phiprime(0,0) = 1; 
    phiprime(1,0) = u1; phiprime(1,1) = p0; 
    phiprime(2,0) = u2; phiprime(2,2) = p0; 
    phiprime(3,0) = u3; phiprime(3,3) = p0; 
    phiprime(4,0) = (u1*u1); phiprime(4,4) = 1; phiprime(4,1) = 2*p0*u1; 
    phiprime(5,0) = u1*u2; phiprime(5,5) = 1; phiprime(5,1) = p0*u2; phiprime(5,2) = p0*u1; 
    phiprime(6,0) = u1*u3; phiprime(6,6) = 1; phiprime(6,1) = p0*u3; phiprime(6,3) = p0*u1; 
    phiprime(7,0) = (u2*u2); phiprime(7,7) = 1; phiprime(7,2) = 2*p0*u2; 
    phiprime(8,0) = u2*u3; phiprime(8,8) = 1; phiprime(8,2) = p0*u3; phiprime(8,3) = p0*u2; 
    phiprime(9,0) = (u3*u3); phiprime(9,9) = 1; phiprime(9,3) = 2*p0*u3; 
    phiprime(10,0) = (u1*u1*u1); phiprime(10,4) = 3*u1; phiprime(10,10) = 1; phiprime(10,1) = 3*p11 + 3*p0*(u1*u1); 
    phiprime(11,0) = (u1*u1)*u2; phiprime(11,4) = u2; phiprime(11,5) = 2*u1; phiprime(11,11) = 1; phiprime(11,1) = 2*p12 + 2*p0*u1*u2; phiprime(11,2) = p11 + p0*(u1*u1); 
    phiprime(12,0) = (u1*u1)*u3; phiprime(12,4) = u3; phiprime(12,6) = 2*u1; phiprime(12,12) = 1; phiprime(12,1) = 2*p13 + 2*p0*u1*u3; phiprime(12,3) = p11 + p0*(u1*u1); 
    phiprime(13,0) = u1*(u2*u2); phiprime(13,5) = 2*u2; phiprime(13,7) = u1; phiprime(13,13) = 1; phiprime(13,1) = p22 + p0*(u2*u2); phiprime(13,2) = 2*p12 + 2*p0*u1*u2; 
    phiprime(14,0) = u1*u2*u3; phiprime(14,5) = u3; phiprime(14,6) = u2; phiprime(14,8) = u1; phiprime(14,14) = 1; phiprime(14,1) = p23 + p0*u2*u3; phiprime(14,2) = p13 + p0*u1*u3; phiprime(14,3) = p12 + p0*u1*u2; 
    phiprime(15,0) = u1*(u3*u3); phiprime(15,6) = 2*u3; phiprime(15,9) = u1; phiprime(15,15) = 1; phiprime(15,1) = p33 + p0*(u3*u3); phiprime(15,3) = 2*p13 + 2*p0*u1*u3; 
    phiprime(16,0) = (u2*u2*u2); phiprime(16,7) = 3*u2; phiprime(16,16) = 1; phiprime(16,2) = 3*p22 + 3*p0*(u2*u2); 
    phiprime(17,0) = (u2*u2)*u3; phiprime(17,7) = u3; phiprime(17,8) = 2*u2; phiprime(17,17) = 1; phiprime(17,2) = 2*p23 + 2*p0*u2*u3; phiprime(17,3) = p22 + p0*(u2*u2); 
    phiprime(18,0) = u2*(u3*u3); phiprime(18,8) = 2*u3; phiprime(18,9) = u2; phiprime(18,18) = 1; phiprime(18,2) = p33 + p0*(u3*u3); phiprime(18,3) = 2*p23 + 2*p0*u2*u3; 
    phiprime(19,0) = (u3*u3*u3); phiprime(19,9) = 3*u3; phiprime(19,19) = 1; phiprime(19,3) = 3*p33 + 3*p0*(u3*u3); 

    Lucee::Matrix<double> reigs(20,20);
    Lucee::Vector<double> evi(20);
    Lucee::Matrix<double> w(20,20);
    Lucee::Vector<double> eigvals(20);

    eigRight(A, eigvals, evi, reigs); 
    w = accumulate(w, phiprime, reigs); 
    solve(reigs, phiJump); // phiJump now contains the left eigenvector projection    
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < 20; j++) {
        waves(j,i) = phiJump(i)*w(j,i);
      }
    }

    // output speeds: 
    for (int i = 0; i < 20; ++i) { 
      s[i] = eigvals(i);
      /*      if (fabs(evi(i)) > 1e-3) {
        std::cout << i <<"'th imaginary part in terms of thermal speed : " << evi(i)/sqrt(p11/p0) << "\n";
        }*/
    }

  }

  void TwentyMomentEquation::wavesLax(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& jump, const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr, Lucee::Matrix<double>& waves, Lucee::FieldPtr<double>& s)
  {
    if (useIntermediateWave)
    {
// this uses the HLLE technique to construct an intermediate state
      std::vector<const double*> auxVars;
      double fl[20], fr[20];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[0]+sr[0]);
      s[1] = 0.5*(sl[1]+sr[1]);

// compute intermediate HLLE state
      double sdiff1 = 1/(s[1]-s[0]);
      double qHHLE[20];
      for (unsigned i=0; i<20; ++i)
        qHHLE[i] = (s[1]*qr[i]-s[0]*ql[i]+fl[i]-fr[i])*sdiff1;

// compute waves
      for (unsigned i=0; i<20; ++i)
      {
        waves(i,0) = qHHLE[i]-ql[i];
        waves(i,1) = qr[i]-qHHLE[i];
      }
    }
    else
    {
// use a single wave
      for (unsigned i=0; i<20; ++i)
        waves(i,0) = qr[i]-ql[i];
      
      double sl[2], sr[2];
      speeds(c, &ql[0], sl);
      speeds(c, &qr[0], sr);
      s[0] = 0.5*(sl[1]+sr[1]);
    }
  }

  void TwentyMomentEquation::qFluctuations(const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& ql, const Lucee::ConstFieldPtr<double>& qr, const Lucee::Matrix<double>& waves, const Lucee::FieldPtr<double>& s, Lucee::FieldPtr<double>& amdq, Lucee::FieldPtr<double>& apdq)
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
      double fl[20], fr[20];
      flux(c, &ql[0], auxVars, fl);
      flux(c, &qr[0], auxVars, fr);

      double absMaxs = std::max(maxAbsSpeed(c, &ql[0]), maxAbsSpeed(c, &qr[0]));
      double amdqL[20], apdqL[20];
      for (unsigned m=0; m<20; ++m)
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

  bool TwentyMomentEquation::isInvariantDomain(const double* q) const
  {
      double rho = q[RHO];
      if(rho <= 0.0){
        std::cout <<"neg density \n";
          return false;
      }

      double u = q[U1] / rho;
      double v = q[U2] / rho;
      double w = q[U3] / rho;

      if((q[P11] - rho * u * u) <= 0) {
        //        std::cout << "neg p11\n";
          return false;
      }
      if((q[P22] - rho * v * v) <= 0){
        //        std::cout << "neg p22 \n";
          return false;
      }
      if((q[P33] - rho * w * w) <= 0) {
        //        std::cout << "neg p33\n";
          return false;
      }
      return true;
  }
}
