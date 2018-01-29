/**
 * @file	LcZonalAverageCalc3D.h
 *
 * @brief	Updater to compute zonal average of a ndim field/distf.
 * Output will be a ndim-1 field.
 */

#ifndef LC_GENERAL_ZONAL_AVERAGE_H
#define LC_GENERAL_ZONAL_AVERAGE_H

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
 * Updater to compute zonal average of a ndim field/distf.
 */
  template <unsigned NDIM>
  class GeneralZonalAverage : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned SNC = Lucee::CDIM<NDIM>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned ANC = Lucee::CDIM<NDIM-1>::N;

/** Create new solver */
      GeneralZonalAverage();
/** Destructor */
      virtual ~GeneralZonalAverage();

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
/** Pointer to nodal source space (NDIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM> *sourceBasis;
/** Pointer to nodal integrated space (ADIM) basis functions to use */
      Lucee::NodalFiniteElementIfc<NDIM-1> *aveBasis;
/** Pointer to nodal source space (NDIM) basis functions to use */
//      Lucee::NodalFiniteElementIfc<NDIM> *targetBasis;
/** Integration direction */
      unsigned intDir;
/** Zeroth moment matrix */
      Eigen::MatrixXd mom0Matrix;
/** Space to store moment data for on a processor */
      Lucee::Field<NDIM-1, double> *momentLocal;
/** Mapping of node in phase-space to node in configuration space */
      std::vector<unsigned> srcAveMap;
/** Boolean for determining if every configuration space node has a phase space node associated with it */
      bool sameAveCoords(unsigned n, unsigned cn, double dxMin,
        unsigned srcDir[NDIM-1], const Eigen::MatrixXd& sourceC, const Eigen::MatrixXd& aveC);
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_GENERAL_ZONAL_AVERAGE_H
