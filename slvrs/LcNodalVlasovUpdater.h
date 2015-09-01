/**
 * @file	LcNodalVlasovUpdater.h
 *
 * @brief	Updater to solve Vlasov equations with nodal DG scheme.
 */

#ifndef LC_NODAL_VLASOV_UPDATER_H
#define LC_NODAL_VLASOV_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcCDIM.h>
#include <LcField.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to solve Vlasov equation using a nodal discontinous
 * Galerkin scheme.
 */
  template <unsigned CDIM, unsigned VDIM>
  class NodalVlasovUpdater : public Lucee::UpdaterIfc
  {
// Number of components for coordinate arrays etc.
      static const unsigned PNC = Lucee::CDIM<CDIM+VDIM>::N;
      static const unsigned CNC = Lucee::CDIM<CDIM>::N;
      
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new nodal DG solver */
      NodalVlasovUpdater();

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
/** Directions to update */
      std::vector<unsigned> updateDims;
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
/** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> phaseConfMap;
/** Species charge */
      double charge;
/** Species mass */
      double mass;

/**
 * Matrix holder: this class is needed as the Matrix class does not
 * have a default constructor.
 */
      struct MatrixHolder
      {
/** Ctor */
          MatrixHolder() : m(1, 1) {}
/** Matrix data */
          Lucee::Matrix<double> m;
      };

/** Stiffness matrices */
      MatrixHolder stiffMatrix[CDIM+VDIM];
/** Lifting matrices for lower surface */
      MatrixHolder lowerLift[CDIM+VDIM];
/** Lifting matrices for upper surface */
      MatrixHolder upperLift[CDIM+VDIM];

/**
 * Structure to store node numbers on edges.
 */
      struct EdgeNodeNums
      {
/** Node numbers */
          std::vector<int> nums;
      };

/** Vector to store lower node numbers */
      EdgeNodeNums lowerNodeNums[CDIM+VDIM];
/** Vector to store upper node numbers */
      EdgeNodeNums upperNodeNums[CDIM+VDIM];

/**
 * Compute physical flux at all nodes in direction 'dir'.
 *
 * @param dir Direction in which flux is needed.
 * @param pc Coordinates of each node in phase-space.
 * @param distf Distribution function at each node.
 * @param EM Electromagnetic field at each configuration space node.
 * @param flux On output, the flux in direction 'dir' at each node.
 */
      void calcFlux(unsigned dir, const Lucee::Matrix<double>& pc,
        const Lucee::ConstFieldPtr<double>& distf, const Lucee::ConstFieldPtr<double>& EM,
        std::vector<double>& flux);

/**
 * Compute numerical flux in direction 'dir'.
 *
 * @param dir Direction in which flux is needed.
 * @param vCoord Velocity coordinate of node.
 * @param distfl Distribution function at left of node.
 * @param distf Distribution function at right of node.
 * @param eml Electromagnetic field at left of node.
 * @param em Electromagnetic field at right of node.
 * @param flux On output, the flux in direction 'dir'.
 * @return Wave speed.
 */
      double numericalFlux(unsigned dir, double vCoord[],
        double distfl, double distf, const double* eml, const double* em,
        double& flux);

/**
 * Compute matrix-vector multiply. Output vector must be
 * pre-allocated. Note that the computation performed is
 *
 * out = m*mat*vec + v*out
 *
 * @param m Factor for accumulation.
 * @param mat Matrix for multiplication.
 * @param vec Vector for multiplication.
 * @param v Factor for accumulation.
 * @param out On output, holds the product.
 */
      void matVec(double m, const Lucee::Matrix<double>& mat,
        const double* vec, double v, double* out);

/**
 * Check if coordinates of phase-space node and configuration-space
 * node are the same.
 */
      bool sameConfigCoords(unsigned n, unsigned cn, double dxMin,
        const Lucee::Matrix<double>& phaseC, const Lucee::Matrix<double>& confC);

      double getSafeVx(int n, const Lucee::Matrix<double>& pc);
      double getSafeVy(int n, const Lucee::Matrix<double>& pc);
      double getSafeVz(int n, const Lucee::Matrix<double>& pc);
  };
}

#endif // LC_NODAL_VLASOV_UPDATER_H
