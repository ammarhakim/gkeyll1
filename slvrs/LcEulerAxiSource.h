/**
 * @file	LcEulerAxiSource.h
 *
 * @brief	Source for computing Lorentz force on a fluid
 */
#ifndef LC_EULER_AXI_SOURCE_H
#define LC_EULER_AXI_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
/**
 * Source for computing Lorentz force on a fluid
 */
  class EulerAxiSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new Lorentz force evaluator.
 */
      EulerAxiSource();

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
/** Gas gamma */
      double gasGamma;
  };
}

#endif // LC_AXI_SOURCE_H
