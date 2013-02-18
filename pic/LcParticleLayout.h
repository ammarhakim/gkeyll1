/**
 * @file	LcParticleLayout.h
 *
 * @brief	Describes layout of a single particle.
 */

#ifndef LC_PARTICLE_LAYOUT_H
#define LC_PARTICLE_LAYOUT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/**
 * This class is used by custom particle classes to define the layout
 * of space and velocity coordinates, and indicate how much total
 * space they need. One can either derive from this class or use it
 * directly with the appropriate temaplte parameters.
 *
 * @tparam XOFF Offset for first spatial coordinate.
 * @tparam VOFF Offset for first velocity coordinate.
 * @tparam NELEMS Total number of data elements in particle.
 */
  template <unsigned XOFF, unsigned VOFF, unsigned NELEMS>
  class ParticleLayout
  {
    public:
/** Offset for first spatial coordinate */
      static const unsigned positionOffset = XOFF;
/** Offset for first velocity coordinate */
      static const unsigned velocityOffset = VOFF;
/** Total number of data elements in particle */
      static const unsigned numElems = NELEMS;      
  };
}

#endif  // LC_PARTICLE_LAYOUT_H
