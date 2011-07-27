/**
 * @file	LcRteHomogeneousSlab.cpp
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 */

// lucee includes
#include <LcArrayIo.h>
#include <LcHdf5Io.h>
#include <LcLinAlgebra.h>
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcMathLib.h>
#include <LcMatrix.h>
#include <LcRteHomogeneousSlab.h>
#include <LcRtePhaseFunction.h>

// gsl includes
#include <gsl/gsl_math.h>

// std includes
#include <iostream>
#include <sstream>
#include <fstream>

namespace Lucee
{
// set class ID for use in registration system
  const char *RteHomogeneousSlab::id = "RteHomogeneousSlab";

  RteHomogeneousSlab::RteHomogeneousSlab()
    : SolverIfc(RteHomogeneousSlab::id), w(1), mu(1), betal(1),
      irradp(&Lucee::FixedVector<2, unsigned>(1)[0],
        &Lucee::FixedVector<2, int>(0)[0]),
      irradm(&Lucee::FixedVector<2, unsigned>(1)[0],
        &Lucee::FixedVector<2, int>(0)[0]),
      dummymu(1)
  {
    isHalfSpace = false;
  }

  RteHomogeneousSlab::~RteHomogeneousSlab()
  {
    delete radiancep;
    delete radiancem;
  }

  void
  RteHomogeneousSlab::readInput(Lucee::LuaTable& tbl)
  {
// read basic parameters
    L = tbl.getNumber("L");
    N = tbl.getNumber("N");
    mu0 = tbl.getNumber("mu0");
    numModes = tbl.getNumber("numModes");
    albedo = tbl.getNumber("albedo");
    flux = tbl.getNumber("flux");
// check if this is a half-space problem
    isHalfSpace = false;
    if (tbl.hasNumber("isHalfSpace"))
    {
      if (tbl.getNumber("isHalfSpace") == 1)
        isHalfSpace = true;
    }

    if (tbl.hasNumVec("tauRadOut"))
    {
// depths at which radiance output is to be written
      tauRadOut = tbl.getNumVec("tauRadOut");
      if (isHalfSpace == false)
      {
        tau0 = tbl.getNumber("tau0");
// ensure that data is not requested deeper tau0
        for (unsigned k=0; k<tauRadOut.size(); ++k)
          if (tauRadOut[k] > tau0)
          {
            Lucee::Except lce(
              "RteHomogeneousSlab::readInput: output can't be requested at depths greater than tau0 (");
            lce << tau0 << ")" << std::endl;
            throw lce;
          }
      }
    }

    if (tbl.hasNumVec("irradOut"))
    {
// get irradiance moments to compute
      std::vector<double> iro = tbl.getNumVec("irradOut");
      for (unsigned i=0; i<iro.size(); ++i)
        irradOut.push_back(iro[i]);
    }

    if (tbl.hasNumVec("tauIrradOut"))
    {
// depths at which irradiance output is to be written
      tauIrradOut = tbl.getNumVec("tauIrradOut");
      if (isHalfSpace == false)
      {
// ensure that data is not requested deeper tau0
        for (unsigned k=0; k<tauIrradOut.size(); ++k)
          if (tauIrradOut[k] > tau0)
          {
            Lucee::Except lce(
              "RteHomogeneousSlab::readInput: output can't be requested at depths greater than tau0 (");
            lce << tau0 << ")" << std::endl;
            throw lce;
          }
      }
    }

    M = 0; // by default no dummy nodes
// check if there are dummy nodes specified
    if (tbl.hasNumVec("dummyNodes"))
    {
      std::vector<double> dn = tbl.getNumVec("dummyNodes");
      M = dn.size();
      dummymu = Lucee::Vector<double>(M);
      for (int i=0; i<M; ++i)
        dummymu[i] = dn[i];
    }

// read in phase function
    Lucee::RtePhaseFunction& pf 
      = tbl.getObjectAsBase<Lucee::RtePhaseFunction>("phaseFunction");
// get expansion coefficients
    betal = pf.getExpCoeffs(L);

// allocate space for weights and ordinates
    w = Lucee::Vector<double>(N);
    mu = Lucee::Vector<double>(N);
// allocate data for storing radiance
    int lo[2] = {0, 0}, up[2];
    up[0] = tauRadOut.size();  up[1] = numModes;
    Lucee::Region<2, int> rgn(lo, up);
    radiancep = new Lucee::Field<2, double>(rgn, N, 0.0);
    radiancem = new Lucee::Field<2, double>(rgn, N, 0.0);
// allocate data for storing irradiance
    if ((tauIrradOut.size() > 0) && (irradOut.size() > 0))
    {
      int start[2] = {0, 0};
      unsigned shape[2];
      shape[0] = tauIrradOut.size(); 
      shape[1] = irradOut.size();
      irradp = Lucee::Array<2, double>(shape, start);
      irradm = Lucee::Array<2, double>(shape, start);
    }

// compute ordinates and wights for use in Gaussian quadrature
    Lucee::gauleg(N, 0, 1, mu, w);
  }

  void
  RteHomogeneousSlab::initialize()
  {
  }

  int
  RteHomogeneousSlab::advance(double t)
  {
// get hold of log stream
    Lucee::Logger& l = Lucee::Logger::get("lucee.console");
    Lucee::LogStream infoStrm = l.getInfoStream();

// allocate various matrices and vectors
    Lucee::Vector<double> nu(N), Nj(N), A(N), B(N), 
      Lp0(N), Lm0(N), Lpt0(N), Lmt0(N), Qp(N), Qm(N);
    Lucee::Matrix<double> phi_p(N, N), phi_m(N, N);

// solve each azimuthal component of RTE
    for (int m=0; m<numModes; ++m)
    {
      infoStrm << "Solving azimuthal mode " << m << " ..." << std::endl;
// compute Qp, Qm
      qbeam(m, Qp, Qm);
// compute eigensystem of RTE
      calc_RTE_eigensystem(m, phi_p, phi_m, nu);
// compute eigensystem normalization
      get_norms(phi_p, phi_m, Nj);
// compute particular solution at top surface
      particular_solution(0, nu, phi_p, phi_m, Nj, Qp, Qm, Lp0, Lm0);
      if (isHalfSpace)
      {
        B = 0.0; // must not allow for growing solutions
// compute coefficients appearing in homogeneous solution
        calc_A_coeffsInf(nu, phi_p, phi_m, Lp0, Lm0, A);
      }
      else
      {
// compute particular solution at bottom surface
        particular_solution(tau0, nu, phi_p, phi_m, Nj, Qp, Qm, Lpt0, Lmt0);
// compute coefficients appearing in homogeneous solution
        calc_AB_coeffs(nu, phi_p, phi_m, Lp0, Lm0, Lpt0, Lmt0, A, B);
      }

// compute irradiances (only m=0)
      if (m==0)
        calc_irradiances(nu, phi_p, phi_m, Nj, Qp, Qm, A, B);

// now compute radiance at requested depths
      Lucee::FieldPtr<double> radp = radiancep->createPtr();
      Lucee::FieldPtr<double> radm = radiancem->createPtr();
      for (unsigned k=0; k<tauRadOut.size(); ++k)
      {
        radiancep->setPtr(radp, k, m);
        radiancem->setPtr(radm, k, m);
// compute particular solution at this depth
        double tau = tauRadOut[k];
        particular_solution(tau, nu, phi_p, phi_m, Nj, Qp, Qm, Lp0, Lm0);
// copy into radiance
        for (int i=0; i<N; ++i)
        {
          radp[i] = Lp0[i];
          radm[i] = Lm0[i];
        }
// now compute total radiance by adding in homogeneous solution
        for (int j=0; j<N; ++j)
        {
          double t1 = A[j]*exp(-tau/nu[j]);
          double t2 = 0;
          if (isHalfSpace == false)
            t2 = B[j]*exp(-(tau0-tau)/nu[j]);
          for (int i=0; i<N; ++i)
          {
            radp[i] += t1*phi_p(i,j) + t2*phi_m(i,j);
            radm[i] += t1*phi_m(i,j) + t2*phi_p(i,j);
          }
        }
      }
    }

    return 0;
  }

  void
  RteHomogeneousSlab::writeToFile(const std::string& baseName, unsigned d)
  {
    if (d == 0)
// don't write anything at start
      return;
// create HDF5 for storing data
    Lucee::Hdf5Io io(0, 0);
    std::string fn = baseName + ".h5";
    Lucee::IoNodeType fNode = io.createFile(fn);

    Lucee::IoNodeType vn;
// write ordinates and weights
    vn = Lucee::writeToFile(io, fNode, "mu", mu);
    io.writeStrAttribute(vn, "description", "ordinates in [0,1]");
    vn = Lucee::writeToFile(io, fNode, "w", w);
    io.writeStrAttribute(vn, "description", "weights in [0,1]");

// save depths at which irradiances are computed
    Lucee::Vector<double> tauIrradOut_v(tauIrradOut.size());
    for (unsigned i=0; i<tauIrradOut_v.getLength(); ++i)
      tauIrradOut_v[i] = tauIrradOut[i];

    vn = Lucee::writeToFile(io, fNode, "tauIrrad", tauIrradOut_v);
    io.writeStrAttribute(vn, "description", 
      "optical depths at which irradiances are output");

// write irradiances to file
    vn = Lucee::writeToFile(io, fNode, "downward_irradiance", irradp);
    io.writeStrAttribute(vn, "units", "W/m^2");

    vn = Lucee::writeToFile(io, fNode, "upward_irradiance", irradm);
    io.writeStrAttribute(vn, "units", "W/m^2");

// arrays to store radiance
    int start[2] = {0, 0};
    unsigned shape[2];
    shape[0] = numModes; shape[1] = N;
    Lucee::Array<2, double> radmi_p(shape, start);
    Lucee::Array<2, double> radmi_m(shape, start);

    Lucee::FieldPtr<double> radp = radiancep->createPtr();
    Lucee::FieldPtr<double> radm = radiancem->createPtr();
// write all radiances
    for (unsigned k=0; k<tauRadOut.size(); ++k)
    {
      std::ostringstream fnp, fnm;
      fnp << "downward_radiance_" << k;
      fnm << "upward_radiance_" << k;
// copy radiance into arrays
      for (int m=0; m<numModes; ++m)
      {
// downward radiance
        radiancep->setPtr(radp, k, m);
// upward radiance
        radiancem->setPtr(radm, k, m);
        for (int i=0; i<N; ++i)
        {
          radmi_p(m,i) = radp[i];
          radmi_m(m,i) = radm[i];
        }
      }
// write radiances
      vn = Lucee::writeToFile(io, fNode, fnp.str(), radmi_p);
      io.writeAttribute(vn, "opticalDepth", tauRadOut[k]);
      io.writeStrAttribute(vn, "units", "W/m^2/sr");

      vn = Lucee::writeToFile(io, fNode, fnm.str(), radmi_m);
      io.writeAttribute(vn, "opticalDepth", tauRadOut[k]);
      io.writeStrAttribute(vn, "units", "W/m^2/sr");
    }
  }

  void
  RteHomogeneousSlab::restoreFromFile(const std::string& baseName)
  {
  }

  void 
  RteHomogeneousSlab::finalize()
  {
  }

  void
  RteHomogeneousSlab::qbeam(int m, Lucee::Vector<double>& Qp, Lucee::Vector<double>& Qm)
  {
    Lucee::Vector<double> PL(N);

    Qp = 0.0;
    Qm = 0.0;
    for (int l=m; l<=L; ++l)
    {
// compute P_l^m(\mu_0)
      double Plmu0 = Lucee::legendre(l, m, mu0);
// compute P_l^m(\mu_i)
      Lucee::legendre(l, m, mu, PL);
      for (int i=0; i<N; i++)
      { // loop over \mu_i
        Qp[i] += Plmu0*PL[i]*betal[l];
        Qm[i] += Plmu0*PL[i]*pow(-1, l-m)*betal[l];
      }
    }
// multiply by \varpi/2*flux
    double fact = 0.5*albedo*flux;
    Qp.scale(fact);
    Qm.scale(fact);
  }

  void
  RteHomogeneousSlab::calc_RTE_eigensystem(int m, Lucee::Matrix<double>& phi_p,
    Lucee::Matrix<double>& phi_m, Lucee::Vector<double>& nu)
  {
// compute matrices E and F
    Lucee::Matrix<double> E(N,N), F(N,N);
    calc_FE(m, F, E);
// compute F*E
    Lucee::Matrix<double> FE(N,N);
    Lucee::accumulate(FE, F, E);
// compute eigenvalues and eignevectors of F*E
    Lucee::Vector<double> evr(N), evi(N); // real and imaginary
    Lucee::Matrix<double> VR(N,N); // right eigenvectors
    Lucee::eigRight(FE, evr, evi, VR);

// compute eigenvalues of RTE
    for (int i=0; i<N; ++i)
      nu[i] = 1.0/sqrt(evr[i]);

// compute eigenvectors of RTE
    Lucee::Vector<double> mu1(N);
    for (int i=0; i<N; ++i)
      mu1[i] = 1/mu[i];

    E.scaleRows(mu1); // M^{-1}*E
// RETHINK THIS: CAN RESCALE X <- nu[j]*M^{1-}*E*X
    for (int j=0; j<N; ++j)
    {
// get j-th right eigenvector
      Lucee::Vector<double> X = VR.getCol(j);

// compute phi_p
      Lucee::Vector<double> phi = phi_p.getCol(j);
// compute phi <- nu[j]*M^{-1}*E*X
      accumulate(0.0, phi, nu[j], E, X);
// do final update to add in M^{-1}*X and
      for (int i=0; i<N; ++i)
        phi[i] = 0.5*(mu1[i]*X[i] + phi[i]);

// compute phi_m
      phi = phi_m.getCol(j);
// compute phi <- nu[j]*M^{-1}*E*X
      accumulate(0.0, phi, nu[j], E, X);
// do final update to add in M^{-1}*X and
      for (int i=0; i<N; ++i)
        phi[i] = 0.5*(mu1[i]*X[i] - phi[i]);
    }
  }

  void
  RteHomogeneousSlab::calc_FE(int m, Lucee::Matrix<double>& F, Lucee::Matrix<double>& E)
  {
// clear existing matrices
    F = 0.0;
    E = 0.0;
    int sign = 1;
    Lucee::Vector<double> PL(N);
    for (int l=m; l<=L; ++l)
    {
      Lucee::legendre(l, m, mu, PL);
// compute \beta_l P_l^m(\mu_i) P_l^m(\mu_j) and put into proper matrix
      if (sign == 1)
        accumulate(E, betal[l], PL, PL);
      else
        accumulate(F, betal[l], PL, PL);
      sign = -1*sign;
    }
// compute \varpi*W*M^{-1}. Note we do not need 1/2 factor as it
// cancels the accumulation done above.
    Lucee::Vector<double> WMinv(N);
    for (int i=0; i<N; ++i)
      WMinv[i] = -albedo*w[i]/mu[i];
// accumulate into F and E matrices
    E.scaleCols(WMinv);
    F.scaleCols(WMinv);

// do final accumulation of M^{-1}
    for (int i=0; i<N; ++i)
    {
      double mu1 = 1/mu[i];
      E(i,i) = mu1 + E(i,i);
      F(i,i) = mu1 + F(i,i);
    }
  }

  void
  RteHomogeneousSlab::get_norms(const Lucee::Matrix<double>& phi_p,
    const Lucee::Matrix<double>& phi_m, Lucee::Vector<double>& Nj)
  {
    for (int i=0; i<N; ++i)
    {
      double sum = 0;
      for (int k=0; k<N; ++k)
      {
        double t1 = phi_p(k,i); // columns of phi_p are eigenvectors
        double t2 = phi_m(k,i); // columns of phi_m are eigenvectors
        sum += (t1*t1 - t2*t2)*mu[k]*w[k];
      }
      Nj[i] = sum;
    }
  }

  void
  RteHomogeneousSlab::particular_solution(double tau, const Lucee::Vector<double>& nu,
    const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& Nj,
    const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
    Lucee::Vector<double>& Lp_p, Lucee::Vector<double>& Lp_m)
  {
    Lp_p = 0.0;
    Lp_m = 0.0;

    Lucee::Vector<double> As(N), Bs(N);
// compute \script{A}_j and \script{B}_j
    scriptAB(tau, nu, phi_p, phi_m, Nj, Qp, Qm, As, Bs);

// compute particular solution
    for (int i=0; i<N; ++i)
    {
      double sum1 = 0.0;
      double sum2 = 0.0;
      for (int j=0; j<N; j++)
      {
        sum1 += (As[j]*phi_p(i,j) + Bs[j]*phi_m(i,j));
        sum2 += (As[j]*phi_m(i,j) + Bs[j]*phi_p(i,j));
      }
      Lp_p[i] = sum1;
      Lp_m[i] = sum2;
    }
  }

  void
  RteHomogeneousSlab::scriptAB(double tau, const Lucee::Vector<double>& nu,
    const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& Nj,
    const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
    Lucee::Vector<double>& As, Lucee::Vector<double>& Bs)
  {
    As = 0.0;
    Bs = 0.0;
    for (int j=0; j<N; ++j)
    {
      double sum1 = 0.0;
      double sum2 = 0.0;
// factor in front of \script{A}
      double t1 = mu0*nu[j]*Cfunc(tau, nu[j], mu0)/Nj[j];
// factor in front of \script{B}
      double t2 = 0.0;
      if (isHalfSpace)
        t2 = mu0*nu[j]*exp(-tau/mu0)*SfuncInf(nu[j], mu0)/Nj[j];
      else
        t2 = mu0*nu[j]*exp(-tau/mu0)*Sfunc(tau0-tau, nu[j], mu0)/Nj[j];

      for (int i=0; i<N; ++i)
      {
        sum1 += w[i]*(Qp[i]*phi_p(i,j) + Qm[i]*phi_m(i,j));
        sum2 += w[i]*(Qp[i]*phi_m(i,j) + Qm[i]*phi_p(i,j));
      }
      As[j] = t1*sum1;
      Bs[j] = t2*sum2;
    }
  }

  double
  RteHomogeneousSlab::Cfunc(double tau, double x, double y)
  {
    if(fabs(x-y) < 1e-14)
      return tau*exp(-tau/y)/(y*y);
    if(fabs(x) < 1e-14)
      return exp(-tau/y)/y;
    if(fabs(y) < 1e-14)
      return exp(-tau/x)/x;
    
    return (exp(-tau/x)-exp(-tau/y))/(x-y);
  }

  double
  RteHomogeneousSlab::Sfunc(double tau, double x, double y)
  {
    if(fabs(x) < 1e-14)
      return 1/y;
    if(fabs(y) < 1e-14)
      return 1/x;
    
    return (1-exp(-tau/x)*exp(-tau/y))/(x+y);
  }

  double
  RteHomogeneousSlab::SfuncInf(double x, double y)
  {
    if(fabs(x) < 1e-14)
      return 1/y;
    if(fabs(y) < 1e-14)
      return 1/x;
    
    return 1/(x+y);
  }

  void
  RteHomogeneousSlab::calc_AB_coeffs(
    const Lucee::Vector<double>& nu,
    const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& Lp0, const Lucee::Vector<double>& Lm0,
    const Lucee::Vector<double>& Lpt0, const Lucee::Vector<double>& Lmt0,
    Lucee::Vector<double>& A, Lucee::Vector<double>& B)
  {
    Lucee::Matrix<double> BLOCKS(2*N, 2*N), RHS(2*N, 1);

    for (int j=0; j<N; ++j)
    {
      double t = exp(-tau0/nu[j]);
      for (int i=0; i<N; ++i)
      {
// compute blocks from top-surface boundary conditions
        BLOCKS(i,j) = phi_p(i,j);
        BLOCKS(i, j+N) = phi_m(i,j)*t;
// compute blocks from bottom-surface boundary conditions
        BLOCKS(i+N,j) = phi_m(i,j)*t;
        BLOCKS(i+N, j+N) = phi_p(i,j);
      }
    }
    for (int i=0; i<N; ++i)
    {
// compute RHS contribution from top boundary contribution
      RHS(i,0) = -Lp0[i];
// compute RHS contribution from bottom boundary contribution
      RHS(i+N,0) = -Lmt0[i];
    }

// invert system to get solution
    Lucee::solve(BLOCKS, RHS);
// copy solution to A and B
    for (int i=0; i<N; ++i)
    {
      A[i] = RHS(i,0);
      B[i] = RHS(i+N,0);
    }
  }

  void
  RteHomogeneousSlab::calc_A_coeffsInf(const Lucee::Vector<double>& nu,
    const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& Lp0, const Lucee::Vector<double>& Lm0,
    Lucee::Vector<double>& A)
  {
    Lucee::Matrix<double> BLOCKS(N, N), RHS(N, 1);

    for (int j=0; j<N; ++j)
    {
      for (int i=0; i<N; ++i)
      {
// compute blocks from top-surface boundary conditions
        BLOCKS(i,j) = phi_p(i,j);
      }
    }
    for (int i=0; i<N; ++i)
    {
// compute RHS contribution from top boundary contribution
      RHS(i,0) = -Lp0[i];
    }

// invert system to get solution
    Lucee::solve(BLOCKS, RHS);
// copy solution to A and B
    for (int i=0; i<N; ++i)
    {
      A[i] = RHS(i,0);
    }
  }

  void
  RteHomogeneousSlab::calc_irradiances(const Lucee::Vector<double>& nu,
    Lucee::Matrix<double>& phi_p, Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& Nj,
    const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
    const Lucee::Vector<double>& A, const Lucee::Vector<double>& B)
  {
// allocate various arrays
    unsigned nirrad = irradOut.size();
    Lucee::Matrix<double> Rp(N, nirrad), Rm(N, nirrad);
    Lucee::Vector<double> PLW(N), PL(N), As(N), Bs(N);

// compute Rp and Rm for each irradiance requested
    for (unsigned mi=0; mi<nirrad; ++mi)
    {
      int mom = irradOut[mi];
// compute PL
      Lucee::legendre(mom, 0, mu, PL);
// compute PL^T*W
      for (int j=0; j<N; ++j) 
        PLW[j] = PL[j]*w[j];

      for (int j=0; j<N; ++j)
      {
        Lucee::Vector<double> phi_p_j = phi_p.getCol(j);
        Rp(j,mi) = PLW.innerProduct(phi_p_j);
        Lucee::Vector<double> phi_m_j = phi_m.getCol(j);
        Rm(j,mi) = PLW.innerProduct(phi_m_j);
      }
    }
// loop over depths
    for (unsigned depth=0; depth<tauIrradOut.size(); ++depth)
    {
      double tau = tauIrradOut[depth];
// compute script{A} and script{B}
      scriptAB(tau, nu, phi_p, phi_m, Nj, Qp, Qm, As, Bs);
// now compute irradiances
      for (unsigned mi=0; mi<irradOut.size(); ++mi)
      {
        int mom = irradOut[mi];
        double sum1 = 0.0, sum2 = 0.0;
        for (int j=0; j<N; j++)
        {
          double t1 = A[j]*exp(-tau/nu[j]) + As[j];
          double t2 = B[j]*exp(-(tau0-tau)/nu[j]) + Bs[j];
          sum1 += t1*Rp(j,mi) + t2*Rm(j,mi);
          sum2 += t1*Rm(j,mi) + t2*Rp(j,mi);
        }
        double direct = Lucee::legendre(mom, 0, mu0)*exp(-tau/mu0)*flux;
        irradp(depth, mi) = M_PI*(sum1+direct);
        irradm(depth, mi) = M_PI*sum2;
      }
    }
  }

  void
  RteHomogeneousSlab::calc_extended_eigensystem(int m,
    const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
    const Lucee::Vector<double>& nu,
    Lucee::Matrix<double>& hphi_p, Lucee::Matrix<double>& hphi_m)
  {
  }
}
