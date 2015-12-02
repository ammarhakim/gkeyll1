/**
 * @file	LcTwentyMomentFluidSource.h
 *
 * @brief       Compute source terms in 20-moment fluid equations.
 */
#ifndef LC_TWENTY_MOMENT_FLUID_SOURCE_H
#define LC_TWENTY_MOMENT_FLUID_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
/**
 * Source for computing source terms in pressure tensor equations.
 */
  class TwentyMomentFluidSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new pressure tensor source evaluator.
 */
      TwentyMomentFluidSource();

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

/**
 * Compute source Jacobian and store it in supplied output
 * matrix. Derived class method should use the getData() method to get
 * data it needs in computing the sources.
 *
 * @param tm Time at which source is requested.
 * @param loc Coordinate at which source is requested.
 * @param jac On output, source jacobian.
 */
      void getSourceJac(double tm, const double loc[3], Lucee::Matrix<double>& jac);

    private:
/** Charge to mass ratio of fluid species */
      double qbym;
/** Flag to indicate if collisions should be included */
      bool hasCollisions;
/** Value of collision frequency (this is fixed for now) */
      double nu;
  };
}

#endif // LC_TWENTY_MOMENT_FLUID_SOURCE_H
