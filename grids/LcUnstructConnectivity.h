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
/**
 * Class to hold unstructured grid connectivity. This class holds
 * connectivity from elements of dimension d to those of dimension
 * dprime.
 */
  class UnstructConnectivity
  {
    public:
/**
 * Create a new connectivity object that connects elements of
 * dimension d to elements of dimension dprime.
 *
 * @param nd Number of elements of dimension d.
 */
      UnstructConnectivity(unsigned nd);

    private:
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

/** No copying allowed */
      UnstructConnectivity(const UnstructConnectivity&);
/** No assignment allowed */
      UnstructConnectivity& operator=(const UnstructConnectivity&);
  };
}

#endif //  LC_UNSTRUCT_GEOMETERY_H
