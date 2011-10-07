/**
 * @file	LcFaceEdgeCurlUpdater.h
 *
 * @brief	Compute curl on rectangular grids.
 */

#ifndef LC_FACE_EDGE_CURL_UPDATER_H
#define LC_FACE_EDGE_CURL_UPDATER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Updater to compute curl of a vector field on a rectangular
 * grid. The output field is assumed to live on cell faces and the
 * input field (to curl) is assumed to live on cell edges.
 */
  template <unsigned NDIM>
  class FaceEdgeCurlUpdater : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new curl updater.
 */
      FaceEdgeCurlUpdater();

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
/** Factor that multiplies curl */
      double alpha;
/** Speed for use in CFL computation */
      double speed;
/** CFL number */
      double cfl;
/** Extra cells to update outside of interior */
      unsigned ghostUpdates[2];
  };
}

#endif // LC_FACE_EDGE_CURL_UPDATER_H
