/**
 * @file	LcFieldFunctionSource.h
 *
 * @brief	Source for computing source from Lua function.
 */
#ifndef LC_FIELD_FUNCTION_SOURCE_H
#define LC_FIELD_FUNCTION_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
/**
 * Source for computing source from Lua function.
 */
  class FieldFunctionSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new function source evaluator.
 */
      FieldFunctionSource();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Compute sources and store them in supplied output vector. The
 * vector 'src' is pre-allocated. Derived class method should use the
 * getData() method to get data it needs in computing the sources.
 *
 * @param tm Time at which source is requested.
 * @param loc Coordinate at which source is requested.
 * @param src On output, source.
 */
      inline void getSource(double tm, const double loc[3], std::vector<double>& src);

    private:
/** Reference to Lua function */
      int fnRef;
/** Input components to send to Lua function */
      std::vector<unsigned> inpComponents;
  };
}

#endif // LC_FIELD_FUNCTION_SOURCE_H
