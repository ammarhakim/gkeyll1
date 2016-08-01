/**
 * @file	LcIonizationSource.h
 *
 * @brief	Source for computing Lorentz force on a fluid
 */
#ifndef LC_IONIZATION_SOURCE_H
#define LC_IONIZATION_SOURCE_H

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
  class IonizationSource : public Lucee::PointSourceIfc
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;

/** 
 * Create new ionization source evaluator.
 */
      IonizationSource();

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
      // Ionization rate constants
      double ionizationConst;
      double ionizationConstDensityElc;
      double ionizationConstDensityIon;
      double ionizationConstEnergy; // energy is not multiplied by
                                    // density of respective species
      // Gas gamma and permeability (used to get temperature from energy)
      double gasGamma;
      double invPermeability;  
      // masses
      double massElc;
      double massIon;
  };
}

#endif // LC_IONIZATION_SOURCE_H
