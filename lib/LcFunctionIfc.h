/**
 * @file	LcFunctionIfc.h
 *
 * @brief	Interface class for function objects.
 */

#ifndef LC_FUNCTION_IFC_H
#define LC_FUNCTION_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * Base class for function-like objects in Lucee
 */
  class FunctionIfc : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/**
 * Create new empty function object. If this default ctor is called
 * then the size of the input and output vector must be specified
 * using setInpSize() and setOutSize() respectively.
 */
      FunctionIfc();

/**
 * Create new function object with specified sizes of the input and
 * output vector.
 */
      FunctionIfc(unsigned ni, unsigned no);

/**
 * Set size of input vector.
 *
 * @param sz size of input vector.
 */
      void setInpSize(unsigned sz) { inpSz = sz; }

/**
 * Set size of output vector.
 *
 * @param sz of output vector.
 */
      void setOutSize(unsigned sz) { outSz = sz; }

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Evaluate the function and return results.
 *
 * @param inp Input values to use.
 * @return result of evaluation.
 */
      virtual std::vector<double> eval(const std::vector<double>& inp) = 0;

    protected:
/**
 * Return size of input vector.
 *
 * @return size of input vector.
 */
      unsigned getInpSize() const { return inpSz; }

/**
 * Return size of output vector.
 *
 * @return size of output vector.
 */
      unsigned getOutSize() const { return outSz; }

    private:
/** Size of input vector */
      unsigned inpSz;
/** Size of output vector */
      unsigned outSz;
  };
}

#endif // LC_FUNCTION_IFC_H
