/**
 * @file	LcRteHomogeneousSlab.h
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 *
 * @version	$Id: LcRteHomogeneousSlab.h 331 2010-03-10 03:55:02Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_RTE_HOMOGENEOUS_SLAB_H
#define LC_RTE_HOMOGENEOUS_SLAB_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcMatrix.h>
#include <LcSolverIfc.h>
#include <LcVector.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to solve the radiative transfer equation in a homogeneous slab.
 */
  class RteHomogeneousSlab : public Lucee::SolverIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create a new RET solver object.
 */
      RteHomogeneousSlab();

/**
 * Destroy solver object.
 */
      virtual ~RteHomogeneousSlab();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Initialize solver, i.e. setup initial conditions. At the end of
 * this call, the solver should be ready for evolving the solution.
 */
      virtual void initialize();

/**
 * Advance the solution to specified time. Solvers that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of solution.
 */
      virtual int advance(double t);

/**
 * Write solver data to file.
 *
 * @param baseName Base name of output files. This should serve as a
 *   prefix for all output files.
 * @param d Dump number.
 */
      virtual void writeToFile(const std::string& baseName, unsigned d);

/**
 * Restore solver data from file. This is called instead of the
 * initialize() method if the simulation is being restarted.
 *
 * @param baseName Base name of input files. This should serves as a
 *   prefix for all input files.
 */
      virtual void restoreFromFile(const std::string& baseName);

/**
 * Finalize solver: free resources, deallocate memory, close files
 * etc.
 */
      virtual void finalize();
    private:
/** Number of expansion coefficients */
      int L;
/** Number of quadrature points in each hemisphere */
      int N;
/** Number of "dummy" quadrature points in each hemisphere */
      int Nd;
/** Number of azimuthal modes */
      int numModes;
/** Cosine of incidence angle */
      double mu0;
/** Downward irradiance is mu0*pi*F */
      double flux;
/** Flag to indicate if this is a half-space problem */
      bool isHalfSpace;
/** Optical depth */
      double tau0;
/** Albedo of single scattering */
      double albedo;
/** Weights for Gaussian quadrature */
      Lucee::Vector<double> w;
/** Ordinates for Gaussian quadrature */
      Lucee::Vector<double> mu;
/** Phase function expasion coefficients */
      Lucee::Vector<double> betal;
/** Irrdiance moments to compute */
      std::vector<int> irradOut;
/** Optical depths at which radiance output is written */
      std::vector<double> tauRadOut;
/** Optical depths at which irradiance output is written */
      std::vector<double> tauIrradOut;
/** Arrays for storing downward irradiances */
      Lucee::Array<2, double> irradp;
/** Arrays for storing upward irradiances */
      Lucee::Array<2, double> irradm;
/** Array to store downward radiance (ntau X numModes) with N components */
      Lucee::Field<2, double> *radiancep;
/** Array to store upward radiance data (ntau X numModes) with N components */
      Lucee::Field<2, double> *radiancem;
/** Number of dummy nodes */
      int M;
/** Set of dummy nodes */
      Lucee::Vector<double> dummymu;

/**
 * Computes Qp and Qm needed to compute particular solutions.
 *
 * @param m Azimuthal mode number.
 * @param Qp on output, beam source in positive hemisphere.
 * @param Qm on output, beam source in negative hemisphere.
 */
      void qbeam(int m, Lucee::Vector<double>& Qp, Lucee::Vector<double>& Qm);

/**
 * Computes the eigensystem of the RTE.
 *
 * @param m Azimuthal mode number.
 * @param phi_p on output, eigenvectors of the RTE corresponding to nu_j
 * @param phi_m on outout, eigenvectors of the RTE  corresponding to -nu_j
 * @param nu on ouput, eigevalues of the RTE
 */
      void calc_RTE_eigensystem(int m, Lucee::Matrix<double>& phi_p,
        Lucee::Matrix<double>& phi_m, Lucee::Vector<double>& nu);

/**
 * Comppute the matrices F and E.
 *
 * @param m azimuthal mode.
 * @param F on ouput, F matrix.
 * @param E on ouput, E matrix.
 */
      void calc_FE(int m, Lucee::Matrix<double>& F, Lucee::Matrix<double>& E);

/**
 * Compute eigesystem normalization.
 *
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Nj on output, normalization coefficients
 */
      void get_norms(const Lucee::Matrix<double>& phi_p,
        const Lucee::Matrix<double>& phi_m, Lucee::Vector<double>& Nj);

/**
 * Compute particular solution at a given depth.
 *
 * @param tau Depth at which particular solution is needed.
 * @param nu Eigenvalues of RTE.
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Nj Normalization coefficients.
 * @param Qp beam source in positive hemisphere.
 * @param Qm beam source in negative hemisphere.
 * @param Lp_p on output, particular solution in positive hemisphere.
 * @param Lp_m on output, particular solution in negative hemisphere.
 */
      void particular_solution(double tau, const Lucee::Vector<double>& nu,
        const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& Nj,
        const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
        Lucee::Vector<double>& Lp_p, Lucee::Vector<double>& Lp_m);

/**
 * Compute functions that appear in the particular solution.
 *
 * @param tau Depth at which functions are to be evaluated.
 * @param nu Eigenvalues of RTE.
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Nj Normalization coefficients.
 * @param Qp beam source in positive hemisphere.
 * @param Qm beam source in negative hemisphere.
 * @param As on output, \script{A}.
 * @param Bs on output, \script{B}.
 */
      void scriptAB(double tau, const Lucee::Vector<double>& nu,
        const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& Nj,
        const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
        Lucee::Vector<double>& As, Lucee::Vector<double>& Bs);

/**
 * C(\tau; x, y) appearing in the particular solution.
 *
 * @param tau Depth.
 * @param x parameter in function.
 * @param y parameter in function.
 */
      double Cfunc(double tau, double x, double y);

/**
 * S(\tau; x, y) appearing in the particular solution.
 *
 * @param tau Depth.
 * @param x parameter in function.
 * @param y parameter in function.
 */
      double Sfunc(double tau, double x, double y);

/**
 * S(x, y) appearing in the particular solution for a half-space
 * problem.
 *
 * @param x parameter in function.
 * @param y parameter in function.
 */
      double SfuncInf(double x, double y);

/**
 * Compute the coefficients appearing in the homogeneous solution.
 *
 * @param nu Eigenvalues of RTE.
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Lp0 Particular solution in [0,1] on top surface.
 * @param Lm0 Particular solution in [0,-1] on top surface.
 * @param Lpt0 Particular solution in [0,1] on bottom surface.
 * @param Lmt0 Particular solution in [0,-1] on bottom surface.
 * @param A on output, A coefficients in homogeneous solution.
 * @param B on output, B coefficients in homogeneous solution.
 */
      void calc_AB_coeffs(const Lucee::Vector<double>& nu,
        const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& Lp0, const Lucee::Vector<double>& Lm0,
        const Lucee::Vector<double>& Lpt0, const Lucee::Vector<double>& Lmt0,
        Lucee::Vector<double>& A, Lucee::Vector<double>& B);

/**
 * Compute the coefficients appearing in the homogeneous solution for
 * a half-space problem.
 *
 * @param nu Eigenvalues of RTE.
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Lp0 Particular solution in [0,1] on top surface.
 * @param Lm0 Particular solution in [0,-1] on top surface.
 * @param A on output, A coefficients in homogeneous solution.
 */
      void calc_A_coeffsInf(const Lucee::Vector<double>& nu,
        const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& Lp0, const Lucee::Vector<double>& Lm0,
        Lucee::Vector<double>& A);

/**
 * Compute required irradiances at specified depths
 *
 * @param nu Eigenvalues of RTE.
 * @param phi_p Eigenvectors correspoding to +ve eigenvalues.
 * @param phi_m Eigenvectors correspoding to -ve eigenvalues.
 * @param Nj Normalization coefficients.
 * @param Qp beam source in positive hemisphere.
 * @param Qm beam source in negative hemisphere.
 * @param A expansion coefficients for homogeneous solution.
 * @param B expansion coefficients for homogeneous solution.
 */
      void calc_irradiances(const Lucee::Vector<double>& nu,
        Lucee::Matrix<double>& phi_p, Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& Nj,
        const Lucee::Vector<double>& Qp, const Lucee::Vector<double>& Qm,
        const Lucee::Vector<double>& A, const Lucee::Vector<double>& B);

/**
 * Calculate the extended eigensystem. Note that only the first N
 * extended eigenvectors are computed.
 *
 * @param m Azimuthal mode number.
 * @param phi_p eigenvectors of the RTE corresponding to nu_j
 * @param phi_m eigenvectors of the RTE  corresponding to -nu_j
 * @param nu eigenvalues of the RTE
 * @param hphi_p on output, extended eigenvectors of the RTE corresponding to nu_j
 * @param hphi_m on outout, extended eigenvectors of the RTE  corresponding to -nu_j
 */
      void calc_extended_eigensystem(int m, 
        const Lucee::Matrix<double>& phi_p, const Lucee::Matrix<double>& phi_m,
        const Lucee::Vector<double>& nu,
        Lucee::Matrix<double>& hphi_p, Lucee::Matrix<double>& hphi_m);
        
  };
}

#endif // LC_RTE_HOMOGENEOUS_SLAB_H
