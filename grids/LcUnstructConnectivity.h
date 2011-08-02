/**
 * @file	LcUnstructConnectivity.h
 *
 * @brief	Class holding connectivity information.
 */

#ifndef LC_UNSTRUCT_CONNECTIVITY_H
#define LC_UNSTRUCT_CONNECTIVITY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

namespace Lucee
{
/**
 * Class to hold unstructured grid connectivity. This class holds
 * connectivity from elements of dimension d to those of dimension
 * dprime.
 */
  class UnstructConnectivity
  {
    public:
/**
 * Create empty connectivity object. The reset method must be called
 * to allocate memory to store the connectivity.
 */
      UnstructConnectivity();

/**
 * Reset connectivity object to connects elements of dimension d to
 * elements of dimension dprime. Calling this method will clear all
 * data in the object.
 *
 * @param nd Number of elements of dimension d.
 */
      void reset(unsigned nd);

/** Indices of grid element dprime */
      std::vector<int> indices;
/** Offsets into 'indices' array: size is numElem(d)+1. The
 * connections of element n of dimension d are stored in indices[j],
 * where offset[n] <= j < offset[n+1].
 */
      std::vector<unsigned> offsets;
  };
}

#endif //  LC_UNSTRUCT_GEOMETERY_H
