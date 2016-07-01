/**
 * @file	LcBoltzmannCollisionOperator.h
 *
 * @brief	Updater to compute a collision operator for RHS of the Bolztmann 
 *     as defined in Meier & Shumlak (2012) [http://dx.doi.org/10.1063/1.4736975]
 */

#ifndef LC_BOLTZMANN_COLLISION_OPERATOR_H
#define LC_BOLTZMANN_COLLISION_OPERATOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>
#include <LcCDIM.h>

// eigen includes
#include <Eigen/LU>

namespace Lucee
{
/**
 * Updater to compute a collision operator for the Boltzmann equation
 */
  template <unsigned CDIM, unsigned VDIM>
  class BoltzmannCollisionOperator : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned CNC = Lucee::CDIM<CDIM>::N;

/** Create new modal DG solver in CDIM dimensions */
      BoltzmannCollisionOperator();

/**  Destructor */
      virtual ~BoltzmannCollisionOperator();
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
/** Pointer to nodal phase space (CDIM+VDIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM+VDIM> *phaseBasis;
/** Pointer to nodal configuration space (CDIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<CDIM> *confBasis;
/** Boolean for determining if every configuration space node has a phase space node associated with it */
      bool sameConfigCoords(unsigned n, unsigned cn, double dxMin,
        const Eigen::MatrixXd& phaseC, const Eigen::MatrixXd& confC);
/** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> phaseConfMap;
/** Moment to compute */
      unsigned calcMom;
/** Space to store moment data for on a processor */
      Lucee::Field<CDIM, double> *momentLocal;
/** Zeroth moment matrix */
      Eigen::MatrixXd mom0Matrix;
/** First moment matrix */
      std::vector<Eigen::MatrixXd> mom1Matrix;
/** Second moment matrix */
      std::vector<Eigen::MatrixXd> mom2Matrix;
/** Number of moments being computed (vector or tensor) */
      unsigned nMom;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_BOLTZMANN_COLLISION_OPERATOR_H
