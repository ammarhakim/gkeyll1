/**
 * @file	LcParticleSpecies.cpp
 *
 * @brief	Container class to store particles
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcParticleSpecies.h>
#include <LcStandardParticle.h>

namespace Lucee
{
// names used in registration system
  template <> const char *ParticleSpecies<float, StandardParticle<float> >::id = "ParticleSpeciesFlt";
  template <> const char *ParticleSpecies<double, StandardParticle<double> >::id = "ParticleSpecies";

  template <typename REAL, class PTCL>
  ParticleSpecies<REAL, PTCL>::ParticleSpecies()
    : DataStructIfc()
  {
  }

  template <typename REAL, class PTCL>
  ParticleSpecies<REAL, PTCL>::~ParticleSpecies()
  {
  }

  template <typename REAL, class PTCL>
  void
  ParticleSpecies<REAL, PTCL>::readInput(Lucee::LuaTable& tbl)
  {
    DataStructIfc::readInput(tbl);
  }

  template <typename REAL, class PTCL>
  TxIoNodeType
  ParticleSpecies<REAL, PTCL>::writeToFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
    return node;
  }

  template <typename REAL, class PTCL>
  TxIoNodeType 
  ParticleSpecies<REAL, PTCL>::readFromFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
    return node;
  }

  template <typename REAL, class PTCL>
  void
  ParticleSpecies<REAL, PTCL>::sync()
  {
  }

// instantiations
  template class ParticleSpecies<float, StandardParticle<float> >;
  template class ParticleSpecies<double, StandardParticle<double> >;
}
