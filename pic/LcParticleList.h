/**
 * @file	LcParticleList.h
 *
 * @brief	Class to store a list of particles.
 */

#ifndef LC_PARTICLE_LIST_H
#define LC_PARTICLE_LIST_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcParticleBase.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * This class allows storing a list of particles, and provides
 * iterators to manipulate them. Inserting particles into the list has
 * asymptotically constant time cost.
 *
 * @tparam PARTICLE Particle class to store in list.
 */
  template <class PARTICLE>
  class ParticleList
  {
    public:
/** Particle list iterator type */
      typedef typename std::vector<PARTICLE>::iterator iterator;
/** Particle list constant-iterator type */
      typedef typename std::vector<PARTICLE>::const_iterator const_iterator;

/**
 * An empty list of particles.
 */
      ParticleList()
      { }

/**
 * Construct a new list of particles. Note that no particles are
 * actually created, but space is reserved for 'np' particles.
 *
 * @param np Reserve space for these many particles.
 */
      ParticleList(unsigned np)
      {
        particleList.reserve(np);
      }

/**
 * Reserve number of particles to store in list.
 *
 * @param nz New size of list.
 */
      void reserve(unsigned nz)
      { 
        particleList.reserve(nz);
      }

/**
 * Number of particles in list.
 *
 * @return Number of particles in list.
 */
      unsigned numParticles() const 
      { 
        return particleList.size(); 
      }

/**
 * Current capacity of the list. This is the nummber of particles that
 * can be inserted before a reallocation is performed.
 *
 * @return Number of particles in list.
 */
      int excessCapacity() const 
      { 
        return particleList.capacity()-numParticles();
      }

/**
 * Add new particle to list.
 *
 * @param ptcl Particle to add.
 */
      void addParticle(const PARTICLE& ptcl)
      {
        particleList.push_back(ptcl);
      }
/**
 * Delete particle at specified index location.
 *
 * @param idx Index location of particle to delete.
 */
      void deleteParticle(unsigned idx)
      {
        unsigned last = particleList.size()-1;
// copy last particle to this location
        particleList[idx].copyParticle(particleList[last]);
// delete last particle
        particleList.pop_back();
      }

/**
 * Delete particle specified by its iterator.
 *
 * @param itr Iterator to particle for delection.
 */
      void deleteParticle(iterator itr)
      {
        unsigned last = particleList.size()-1;
// copy last particle to this location
        itr->copyParticle(particleList[last]);
// delete last particle
        particleList.pop_back();
      }

/**
 * Starting iterator to particle in list.
 *
 * @return Starting iterator in particle list.
 */
      iterator begin() 
      { 
        return particleList.begin(); 
      }

/**
 * Starting constant iterator to particle in list.
 *
 * @return Starting iterator in particle list.
 */
      const_iterator begin() const 
      { 
        return particleList.begin(); 
      }

/**
 * End iterator to particle in list.
 *
 * @return Starting iterator in particle list.
 */
      iterator end() 
      { 
        return particleList.end(); 
      }

/**
 * End constant iterator to particle in list.
 *
 * @return Starting iterator in particle list.
 */
      const_iterator end() const 
      { 
        return particleList.end(); 
      }

/**
 * Return reference to specified particle.
 *
 * @param loc Location of particle in list.
 * @param reference to particle.
 */
      PARTICLE& operator[](unsigned loc)
      { 
        return particleList[loc]; 
      }

/**
 * Return constant reference to specified particle.
 *
 * @param loc Location of particle in list.
 * @param reference to particle.
 */
      const PARTICLE& operator[](unsigned loc) const
      { 
        return particleList[loc]; 
      }

    private:
/** List of particles */
      std::vector<PARTICLE> particleList;
  };
}

#endif // LC_PARTICLE_LIST_H
