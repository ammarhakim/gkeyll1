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
// \sum_{l=m}^L \beta_l P_l^m(\mu_0) P_l^m(\mu)
    for (unsigned l=m; l<=L; ++l)
    {
      double Plmu0 = Lucee::legendre(l, m, mu0);
      Lucee::legendre(l, m, mu, PL);
      for (unsigned i=0; i<N; i++)
      {
        Qp[i] += PL[i]*betal[l]*Plmu0;
        Qm[i] += PL[i]*pow(-1, l-m)*betal[l]*Plmu0;
      }
    }
    Qp.scale(fact);
    Qm.scale(fact);
  }
}
