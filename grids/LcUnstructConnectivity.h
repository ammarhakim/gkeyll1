/**
 * @file	LcUnstructConnectivity.h
 *
 * @brief	Class holding connectivity information.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
// forward declare creator class
  template <typename REAL> class UnstructGridCreator;
// forward declare grid class
  template <typename REAL> class UnstructGrid;

/**
 * Class to hold unstructured grid connectivity. This class holds
 * connectivity from elements of dimension d to those of dimension
 * dprime. This class is private and hence can not be accessed
 * directly.
 */
  class UnstructConnectivity
  {
    public:
// declare friends so they can touch our privates
      template <typename REAL> friend class UnstructGridCreator;
      template <typename REAL> friend class UnstructGrid;

/**
 * Create empty connectivity object. The reset method must be called
 * to allocate memory to store the connectivity.
 */
      UnstructConnectivity();

    private:
/**
 * Reset connectivity object to connects elements of dimension d to
 * elements of dimension dprime. Calling this method will clear all
 * data in the object.
 *
 * @param nd Number of elements of dimension d.
 */
      void reset(unsigned nd);

/** Source element dimension */
      unsigned d;
/** Target element dimension */
      unsigned dprime;
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
