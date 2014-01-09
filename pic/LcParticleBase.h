/**
 * @file	LcParticleBase.h
 *
 * @brief	Base class for particles.
 */

#ifndef LC_PARTICLE_BASE_H
#define LC_PARTICLE_BASE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcParticleLayout.h>

// std includes
#include <string.h>

namespace Lucee
{
/**
 * This class holds the basic information for a particle. Derived
 * classes can provide specific methods and data particular to that
 * particle type. Often, it is enough to simply inherit from this
 * class without providing any additional methods in the derived
 * class.
 *
 * @tparam REAL Real type (float, double, ..) for storing data.
 * @tparam LAYOUT ParticleLayout class that describes particle data layout.
 */
  template <typename REAL, class LAYOUT>
  class ParticleBase 
  {
    public:
/** Type of real numbers to use */
      typedef REAL value_type;
/** Number of elements per particle */
      static const unsigned numElems = LAYOUT::numElems;

/**
 * Fetch coordinate data.
 *
 * @param d Direction.
 * @return coordinate in specified direction.
 */
      REAL x(unsigned d) const
      {
        return getLoc(LAYOUT::positionOffset+d);
      }

/**
 * Fetch velocity data.
 *
 * @param d Direction.
 * @return velocity in specified direction.
 */
      REAL v(unsigned d) const
      { 
        return getLoc(LAYOUT::velocityOffset+d);
      }

/**
 * Set coordinate data.
 *
 * @param d Direction.
 * @param xn Value to set.
 */
      void setx(unsigned d, REAL xn)
      { 
        setLoc(LAYOUT::positionOffset+d, xn);
      }

/**
 * Set velcity data.
 *
 * @param d Direction.
 * @param xn Value to set.
 */
      void setv(unsigned d, REAL vn)
      { 
        setLoc(LAYOUT::velocityOffset+d, vn);
      }

/**
 * Increment coordinate.
 *
 * @param d Direction.
 * @param xn Value to increment current coordinate by.
 */
      void incrx(unsigned d, REAL xn)
      { 
        incrLoc(LAYOUT::positionOffset+d, xn);
      }

/**
 * Increment velocity
 *
 * @param d Direction.
 * @param vn Value to increment current velocity by.
 */
      void incrv(unsigned d, REAL vn)
      { 
        incrLoc(LAYOUT::velocityOffset+d, vn);
      }

/**
 * Copy supplied particle into this one.
 *
 * @param srcPtcl Particle data to copy over.
 */
      void copyParticle(const ParticleBase<REAL, LAYOUT>& srcPtcl)
      {
        memcpy(data, srcPtcl.data, numElems*sizeof(REAL));
      }

    protected:
/**
 * Set specified location in particle data array.
 *
 * @param loc Location to set.
 * @param val Value to set.
 */
      void setLoc(unsigned loc, REAL val)
      {
        data[loc] = val;
      }

/**
 * Get specified location in particle data array.
 *
 * @param loc Location to set.
 * @return value at location
 */
      REAL getLoc(unsigned loc) const
      {
        return data[loc];
      }

/**
 * Increment value at speicified location.
 *
 * @param loc Location.
 * @param val value to increment by.
 */
      REAL incrLoc(unsigned loc, REAL val)
      {
        data[loc] += val;
      }

    private:
/** Particle data */
      REAL data[LAYOUT::numElems];
  };
}

#endif // LC_PARTICLE_BASE_H
