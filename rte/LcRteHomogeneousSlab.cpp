/**
 * @file	LcRteHomogeneousSlab.cpp
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 *
 * @version	$Id: LcRteHomogeneousSlab.cpp 331 2010-03-10 03:55:02Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLinAlgebra.h>
#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcMathLib.h>
#include <LcMatrix.h>
#include <LcObjCreator.h>
#include <LcRteHomogeneousSlab.h>
#include <LcRtePhaseFunction.h>

// std includes
#include <iostream>

namespace Lucee
{
// set class ID for use in registration system
  const char *RteHomogeneousSlab::id = "RteHomogeneousSlab";

  RteHomogeneousSlab::RteHomogeneousSlab()
    : SolverIfc(RteHomogeneousSlab::id), w(1), mu(1), betal(1)
  {
  }

  RteHomogeneousSlab::~RteHomogeneousSlab()
  {
  }

  void
  RteHomogeneousSlab::readInput(Lucee::LuaTable& tbl)
  {
// read basic parameters
    L = tbl.getNumber("L");
    N = tbl.getNumber("N");
    mu0 = tbl.getNumber("mu0");
    tau0 = tbl.getNumber("tau0");
    numModes = tbl.getNumber("numModes");
    albedo = tbl.getNumber("albedo");
    flux = tbl.getNumber("flux");

// read in phase function
    Lucee::LuaTable pfTbl = tbl.getTable("phaseFunction");
    std::string pfKind = pfTbl.getKind(); // kind of phase function
    Lucee::RtePhaseFunction *pf =
      Lucee::ObjCreator<Lucee::RtePhaseFunction>::getNew(pfKind);
// initialize PF from its table
    pf->readInput(pfTbl);
// get expansion coefficients
    betal = pf->getExpCoeffs(L);
    delete pf; // no need for PF object
  }

  void 
  RteHomogeneousSlab::buildData()
  {
// allocate space for weights and ordinates
    w = Lucee::Vector<double>(N);
    mu = Lucee::Vector<double>(N);
  }

  void 
  RteHomogeneousSlab::buildAlgorithms()
  {
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
      Lp(N), Lm(N), Qp(N), Qm(N);
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

// assemble BLOCKS and RHS matrices to determine unknown expansion coefficients
      Lucee::Matrix<double> BLOCK(2*N, 2*N);
      Lucee::Matrix<double> RHS(2*N, 1);

// compute particular solution at top surface
      particular_solution(0, nu, phi_p, phi_m, Nj, Qp, Qm, Lp, Lm);

// compute contribution to BLOCKS and RHS from top-surface BCs
      Lucee::Matrix<double> b1 = BLOCK.getView(0, N, 0, N);
      Lucee::Matrix<double> b2 = BLOCK.getView(0, N, N, 2*N);
      Lucee::Matrix<double> rhsb = RHS.getView(0, N, 0, 1);

// compute particular solution at bottom surface
      particular_solution(tau0, nu, phi_p, phi_m, Nj, Qp, Qm, Lp, Lm);

// compute contribution to BLOCKS and RHS from bottom-surface BCs
      b1 = BLOCK.getView(N, 2*N, 0, N);
      b2 = BLOCK.getView(N, 2*N, N, 2*N);
      rhsb = RHS.getView(N, 2*N, 0, 1);
    }

    return 0;
  }

  void
  RteHomogeneousSlab::writeToFile(const std::string& baseName, unsigned d) const
  {
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
        Qp[i] += PL[i]*betal[l];
        Qm[i] += PL[i]*pow(-1, l-m)*betal[l];
      }
      Qp.scale(Plmu0);
      Qm.scale(Plmu0);
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
    Lucee::Matrix<double> VL(N,N), VR(N,N); // left and right
    Lucee::eig(FE, evr, evi, VL, VL);

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
    for (int l=m; l<=L; ++l)
    {
      Lucee::Vector<double> PL(N);
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
      E(i,i) += mu1;
      F(i,i) += mu1;
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
    double tau1 = tau0-tau;
    for (int j=0; j<N; ++j)
    {
      double sum1 = 0.0;
      double sum2 = 0.0;
// factor in front of \script{A}
      double t1 = mu0*nu[j]*flux*Cfunc(tau, nu[j], mu0)/Nj[j];
// factor in front of \script{B}
      double t2 = mu0*nu[j]*flux*exp(-tau/mu0)*Sfunc(tau1, nu[j], mu0)/Nj[j];

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
}
