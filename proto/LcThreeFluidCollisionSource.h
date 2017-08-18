/**
 * @file	LcThreeFluidCollisionSource.h
 *
 * @brief	Source for three fluid collisions [Meier & Shumlak, 2012]
 */
#ifndef LC_THREE_FLUID_COLLISION_SOURCE_H
#define LC_THREE_FLUID_COLLISION_SOURCE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPointSourceIfc.h>

namespace Lucee
{
/**
 * Source for computing three fluid (electrons, ions, and neutrals)
 * collisions [Meier & Shumlak, Phys. Plasmas, 2012]
 */
  class ThreeFluidCollisionSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new ionization source evaluator.
 */
      ThreeFluidCollisionSource();

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
      // Voronov ionization fit parameters
      double phi;
      double voronovA;
      double voronovP;
      double voronovX;
      double voronovK;

      double gasGamma;
      double massElc, massIon, massNeut;
      double chargeElc, chargeIon;
  };
}

#endif // LC_THREE_FLUID_COLLISION_SOURCE_H
