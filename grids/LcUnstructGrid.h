/**
 * @file	LcUnstructGrid.h
 *
 * @brief	Unstructured grid class.
 *
 * @version	$Id$
 */

#ifndef LC_UNSTRUCT_GRID_H
#define LC_UNSTRUCT_GRID_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructConnectivity.h>
#include <LcUnstructGeometry.h>

// std includes
#include <map>
#include <vector>

namespace Lucee
{
/**
 * Class to represent an unstructured grid.
 */
  template <typename REAL>
  class UnstructGrid
  {
    public:

    private:
/** Dimension of grid */
      unsigned ndim;
/** Geometry information (assume 3D) */
      Lucee::UnstructGeometry<3, REAL> geometry;
/** Location d+3*dprime indicates where in 'connectivity' d->dprime connection is stored */
      std::vector<unsigned> ddprime;
/** Connectivity information (potentally has 9 elements) */
      std::vector<Lucee::UnstructConnectivity> connectivity;
  };
}

#endif // LC_UNSTRUCT_GRID_H
