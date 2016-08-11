/**
 * @file	LcParticleSpecies.h
 *
 * @brief	Container class to store particles
 */

#ifndef LC_PARTICLE_SPECIES_H
#define LC_PARTICLE_SPECIES_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>

namespace Lucee
{
/**
 * Container class for a species of particles, living on a grid.
 *
 * @tparam REAL Real type (float, double, ..) for storing data.
 * @tparam PTCL Particle class, providing methods needed in accessing particle data
 */
  template <typename REAL, class PTCL>
  class ParticleSpecies : public Lucee::DataStructIfc
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/** Construct a new particle species object */
      ParticleSpecies();

/** Destroy dataStruct */
      virtual ~ParticleSpecies();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      void readInput(Lucee::LuaTable& tbl);

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the data-struct as it should appear in output.
 * @return node to which data was written.
 */
      TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Read dataStruct from given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to read data from.
 * @param nm Name of the data-struct as it appears in input.
 * @return node from which data was read.
 */
      TxIoNodeType readFromFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Synchronise data in ghost cells, copying skin data into neighbor
 * ghost cells.
 */
      void sync();

    private:
  };
}

#endif // LC_PARTICLE_SPECIES_H
