/**
 * @file	LcLinCombiner.h
 *
 * @brief	Linear combiner updater.
 */

#ifndef LC_LIN_COMBINER_H
#define LC_LIN_COMBINER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUpdaterIfc.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Updater to perform linear combination of input fields.
 */
  template <unsigned NDIM>
  class LinCombiner : public Lucee::UpdaterIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/**
 * Create new combiner.
 */
      LinCombiner();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

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

/**
 * Set coefficients for linear combination.
 *
 * @param coeff Coefficients for linear combination.
 */
      void setCoeff(const std::vector<double>& coeff);

    private:
/** Coefficients for linear combination */
      std::vector<double> coeff;
  };
}

#endif // LC_LIN_COMBINER_H
