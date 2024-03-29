/**
 * @file	LcElectrostaticContPhiUpdater.h
 *
 * @brief	Updater to compute phi using a fixed value of k_perp*rho_s
 */

#ifndef LC_ELECTROSTATIC_CONT_PHI_UPDATER_H
#define LC_ELECTROSTATIC_CONT_PHI_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// petsc includes
#include <petsc.h>

// std includes
#include <vector>

// eigen includes
#include <Eigen/Dense>
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve hyperbolic equations using a nodal discontinous
 * Galerkin scheme.
 */
  class ElectrostaticContPhiUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      ElectrostaticContPhiUpdater();

/** Destroy updater */
      ~ElectrostaticContPhiUpdater();

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
 * Advance the solution to specified time. Updaters that do not have a
 * concept of time should ignore the time parameter.
 *
 * @param t Time to advance the solution to.
 * @return Status of updater.
 */
      Lucee::UpdaterStatus update(double t);

/**
 * Declare the types of input and output variables accepted by this
 * updater. This must be provided by the derived classes for the
 * type-checking to pass. Inside the implementation of this method the
 * derived class must make a sequence of appendInpVarType() and
 * appendOutVarType() calls to declare the input/output data structure
 * types.
 */
      void declareTypes();

    private:
/** Pointer to nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<1> *nodalBasis;
/** Value of k_perp0*rho_s */
      double kPerpTimesRho;
/** Fixed electron temp parameter (eV) */
      double Te0;
/** Flag to determine if cutoff velocities will be used to set phi_s */
      bool useCutoffVelocities;
/** Petsc matrices to store linear operator */
      Mat stiffMat;
/** Petsc vectors for source and initial guess */
      Vec globalSrc, initGuess;
/** Krylov subspace method context */
      KSP ksp;
/** Index set for making copies of solution */
      IS is;
/** Petsc vectors for local copy of data */
      Vec localData;
/** Scatter object to get data onto local processor */
      VecScatter vecSctr;
/** Map of rows to Dirichlet BC values */
      std::map<int, double> rowBcValues;
/** Flags to indicated periodic directions */
      bool periodicFlgs[1];
/** Should source be adjusted? */
      bool adjustSource;
/** Flag to indicate if all directions are periodic */
      bool allPeriodic;
/**
 * Matrix of gaussian quadrature locations.
 * There are three columns by default for (x,y,z)
 * and each row is a different quadrature point for doing integrals.
 */
      Eigen::MatrixXd gaussOrdinates;
/** Weights for gaussian quadrature points */
      std::vector<double> gaussWeights;
/** Vector containing various triple-product basis integrals */
      std::vector<Eigen::MatrixXd> tripleProducts;
/** 
 * Interpolation matrix for bringing quantities from nodal locations to
 * gaussian quadrature points.
 */
      Eigen::MatrixXd interpMatrix;
/** Mass matrix */
      Eigen::MatrixXd massMatrix;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
/**
 * Copy data from Petsc field to Gkeyll fields.
 *
 * @param ptFld Input Petsc field.
 * @param gkFld Output Gkeyll field.
 */
      void copyFromPetscField(Vec ptFld, Lucee::Field<1, double>& gkFld);

/**
 * Copy data from Gkeyll to Petsc fields.
 *
 * @param gkFld Input Gkeyll field.
 * @param ptFld Output Petsc field.
 */
      void copyFromGkeyllField(const Lucee::Field<1, double>& gkFld, Vec ptFld);

      void DoPetscAssembly(Mat& mat, bool isFinal = true);

/**
 ** Compute integral of field over the complete domain.
 **
 ** @param fld Field to integrate.
 ** @param shareFlag Flag to indicate if nodes are shared.
 ** @return integral of field over domain.
 **/
      double getFieldIntegral(const Lucee::Field<1, double>& fld, bool shareFlag);
  };
}

#endif // LC_ELECTROSTATIC_CONT_PHI_UPDATER_H
