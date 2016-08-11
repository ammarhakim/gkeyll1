/**
 * @file	LcEigenNodalVlasovUpdater.h
 *
 * @brief	Updater to solve Vlasov equations with nodal DG scheme using Eigen library.
 */

#ifndef LC_EIGEN_NODAL_VLASOV_UPDATER_H
#define LC_EIGEN_NODAL_VLASOV_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to solve Vlasov equation using a nodal discontinous
 * Galerkin scheme.
 */
  template <unsigned CDIM, unsigned VDIM>
  class EigenNodalVlasovUpdater : public Lucee::UpdaterIfc
  {
// Number of components for coordinate arrays etc.
      static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
      static const unsigned CNC = Lucee::CDIM<CDIM>::N;
      
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      EigenNodalVlasovUpdater();

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
/** Pointer to phase-space basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
/** Pointer to configuration-space basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;
/** CFL number to use */
      double cfl;
/** Maximum CFL number allowed */
      double cflm;
/** Flag to indicate if to only compute increments */
      bool onlyIncrement;
/** Species charge */
      double charge;
/** Species mass */
      double mass;
/** Flag to indicate if one should skip the velocity space sweeps */
      bool skipVelocitySweep;
/** Flag to indicate if we should appky zero-flux BCs in velocity space */
      bool applyZeroFluxBc;
/** Offsets for zero-flux directions along lower edges */
      int lowerZeroFluxOffset[CDIM+VDIM];
/** Offsets for zero-flux directions  along upper edges */
      int upperZeroFluxOffset[CDIM+VDIM];

/**
 * Struct to hold data for Guassian quadrature.
 */
      struct GaussQuadData
      {
/**
 * Reset object.
 * 
 * @param numQuad Numer of quadrature points.
 * @param nlocal Total number of local nodes.
 */
          void reset(int numQuad, int nlocal, int numCoords)
          {
            // allocate memory for various matrices
            weights = Eigen::VectorXd(numQuad);
            interpMat = Eigen::MatrixXd(numQuad, nlocal);

          }
/** Vector of weights */
          Eigen::VectorXd weights;
/** Interpolation matrix */
          Eigen::MatrixXd interpMat;
      };
/** Data for volume quadrature */
      GaussQuadData volQuad;
/** Data for quadrature on each lower face */
      GaussQuadData surfLowerQuad[CDIM+VDIM];
/** Data for quadrature on each upper face */
      GaussQuadData surfUpperQuad[CDIM+VDIM];
/** Inverse of mass matrix */
      Eigen::MatrixXd massMatrixInv;
/** These matrices store the gradients of basis functions evaluated at quadrature points */
      std::vector<Eigen::MatrixXd> bigStoredUpperSurfMatrices;
      std::vector<Eigen::MatrixXd> bigStoredLowerSurfMatrices;
      std::vector<Eigen::MatrixXd> bigStoredVolMatrices;

/**
 * Compute physical flux at all nodes in direction 'dir'.
 *
 * @param dir Direction in which flux is needed.
 * @param pc Coordinates of each node in phase-space.
 * @param distf Distribution function at each node.
 * @param EM Electromagnetic field at each configuration space node.
 * @param flux On output, the flux in direction 'dir' at each node.
 */
      void calcFlux(unsigned dir, const Lucee::Matrix<double>& pc, const Lucee::ConstFieldPtr<double>& EM, Eigen::VectorXd& flux);

/**
 * @param alphaLeft: characteristic velocities at left edge quad points (dot n)
 * @param alphaRight: characteristic velocities at right edge quad points (dot n)
 * @param leftValsAtQuad: left cell values at quad points
 * @param rightValsAtQuad: right cell values at quad points
 * @param numericalFluxAtQuad: resulting numerical flux calculation at quad points
 */
      void computeNumericalFlux(const Eigen::VectorXd& alphaLeft, const Eigen::VectorXd& alphaRight, const Eigen::VectorXd& leftValsAtQuad,
          const Eigen::VectorXd& rightValsAtQuad, Eigen::VectorXd& numericalFluxAtQuad);

/** Boolean for determining if every configuration space node has a phase space node associated with it */
      bool sameConfigCoords(unsigned n, unsigned cn, double dxMin,
        const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC);
/** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> phaseConfMap;

/**
 * Get velocity coordinates "safely", i.e. if there are not enough
 * velocity components, these functions return 0.0.
 */      
      double getSafeVx(int n, const Lucee::Matrix<double>& pc);
      double getSafeVy(int n, const Lucee::Matrix<double>& pc);
      double getSafeVz(int n, const Lucee::Matrix<double>& pc);
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);

// helper to index EM fields at nodes      
      unsigned emidx(unsigned n, unsigned i);      
  };
}

#endif // LC_EIGEN_NODAL_VLASOV_UPDATER_H
