/**
 * @file	LcNeutralDragForceSource.h
 *
 * @brief	Source for computing neutral drag force on a fluid
 */
#ifndef LC_NEUTRAL_DRAG_FORCE_SOURCE_H
#define LC_NEUTRAL_DRAG_FORCE_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
/**
 * Source for computing neutral drag force on a fluid
 */
  class NeutralDragForceSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new neutral drag force evaluator.
 */
      NeutralDragForceSource();

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
      // vector of neutral velocities and collision frequency
      std::vector<double> velocityNeut;
      double nu;
  };
}

#endif // LC_NEUTRAL_DRAG_FORCE_SOURCE_H
