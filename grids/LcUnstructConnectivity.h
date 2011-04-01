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
    private:
/** Source element dimension */
      unsigned d;
/** Target element dimension */
      unsigned dprime;
/** Indices of grid element dprime */
      std::vector<int> indices;
/** Offsets into 'indices' array: size is numElem(d)+1 */
      std::vector<unsigned> offsets;
  };
}

#endif //  LC_UNSTRUCT_GEOMETERY_H
