/**
 * @file	LcRunningAverageOfFieldCalc.h
 *
 * @brief	Using the past N inputs, computes average of field at each node
 */

#ifndef LC_RUNNING_AVERAGE_OF_FIELD_CALC_H
#define LC_RUNNING_AVERAGE_OF_FIELD_CALC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
//#include <LcDynVector.h>
#include <LcField.h>
#include <LcUpdaterIfc.h>

// std includes
#include <list>

namespace Lucee
{
/**
 * Updater to record field in a cell
 */
  template <unsigned NDIM>
  class RunningAverageOfFieldCalc : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** Create new modal DG solver in 1D */
      RunningAverageOfFieldCalc();

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
/** Size of history to keep at every node in local region */
      int sampleSize;
/** Local record */
      std::vector<std::list<std::vector<double> > > localRecord;
  };
}

#endif // LC_RUNNING_AVERAGE_OF_FIELD_CALC_H
