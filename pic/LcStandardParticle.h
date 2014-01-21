/**
 * @file	LcStandardParticle.h
 *
 * @brief	A particle with 3 spatial and 3 velocity coordinates.
 */

#ifndef LC_STANDARD_PARTICLE_H
#define LC_STANDARD_PARTICLE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcParticleBase.h>
#include <LcParticleLayout.h>

namespace Lucee
{
/**
 * Particle class with 3 spatial and 3 velocity coordinates.
 */
  template <typename REAL>
  class StandardParticle : public Lucee::ParticleBase<REAL,
    Lucee::ParticleLayout<0, 3, 6> >
  {
  };
}

#endif // LC_STANDARD_PARTICLE_H
