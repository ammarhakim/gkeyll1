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
    for (unsigned m=0; m<numModes; ++m)
    {
      infoStrm << "Solving azimuthal mode " << m << " ..." << std::endl;
// compute Qp, Qm
      qbeam(m, Qp, Qm);
// compute eigensystem of RTE
      calc_RTE_eigensystem(m, phi_p, phi_m, nu);
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
    double fact = 0.5*albedo;

    Qp = 0.0;
    Qm = 0.0;
    for (unsigned l=m; l<=L; ++l)
    {
// compute P_l^m(\mu_0)
      double Plmu0 = Lucee::legendre(l, m, mu0);
// compute P_l^m(\mu_i)
      Lucee::legendre(l, m, mu, PL);
      for (unsigned i=0; i<N; i++)
      { // loop over \mu_i
        Qp[i] += PL[i]*betal[l];
        Qm[i] += PL[i]*pow(-1, l-m)*betal[l];
      }
      Qp.scale(Plmu0);
      Qm.scale(Plmu0);
    }
// multiply by \varpi/2
    Qp.scale(fact);
    Qm.scale(fact);
  }

  void
  RteHomogeneousSlab::calc_RTE_eigensystem(int m, Lucee::Matrix<double>& phi_p,
    Lucee::Matrix<double>& phi_m, Lucee::Vector<double>& nu)
  {
// compute matrices E and F
    Lucee::Matrix<double> E(N,N), F(N,N);
    calc_FE(F, E);
// compute F*E
    Lucee::Matrix<double> FE(N,N);
    Lucee::accumulate(FE, F, E);
// compute eigenvalues and eignevectors of F*E
    Lucee::Vector<double> evr(N), evi(N); // real and imaginary
    Lucee::Matrix<double> VL(N,N), VR(N,N); // left and right
    Lucee::eig(FE, evr, evi, VL, VL);

// compute eigenvalues of RTE
    for (unsigned i=0; i<N; ++i)
      nu[i] = 1.0/sqrt(evr[i]);

// compute eigenvectors of RTE
  }

  void
  RteHomogeneousSlab::calc_FE(Lucee::Matrix<double>& F, Lucee::Matrix<double>& E)
  {
  }
}
