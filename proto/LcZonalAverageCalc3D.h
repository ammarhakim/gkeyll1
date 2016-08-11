/**
 * @file	LcZonalAverageCalc3D.h
 *
 * @brief	Updater to compute zonal (y) average of a 3d potential.
 * Output will be a 2d field.
 */

#ifndef LC_ZONAL_AVERAGE_CALC_3D_H
#define LC_ZONAL_AVERAGE_CALC_3D_H

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
 * Updater to compute zonal average of phi(x,y,z)
 */
  class ZonalAverageCalc3D : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC3 = Lucee::CDIM<3>::N;
/** Number of components for coordinate arrays etc. */
      static const unsigned NC2 = Lucee::CDIM<2>::N;

/** Create new solver */
      ZonalAverageCalc3D();
/** Destructor */
      virtual ~ZonalAverageCalc3D();

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
/** Pointer to 3D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<3> *nodalBasis3d;
/** Pointer to 2D nodal basis functions to use */
      Lucee::NodalFiniteElementIfc<2> *nodalBasis2d;
/** Zeroth moment matrix */
      Eigen::MatrixXd mom0Matrix;
/** Space to store moment data for on a processor */
      Lucee::Field<2, double> *moment;
/**
 * Copy a Lucee-type matrix to an Eigen-type matrix.
 * No checks are performed to make sure source and destination matrices are
 * of the same size.
 */
      void copyLuceeToEigen(const Lucee::Matrix<double>& sourceMatrix, Eigen::MatrixXd& destinationMatrix);
  };
}

#endif // LC_ZONAL_AVERAGE_CALC_3D_H
