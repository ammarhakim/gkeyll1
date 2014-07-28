/**
 * @file	LcCDIM.h
 *
 * @brief	Classes to provide compile-time size for storing coordinates
 */

#ifndef LC_CDIM_H
#define LC_CDIM_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/**
 * These set of classes give a means of determining number of elements
 * of coordinate-based objects, like physical coordinates, normal and
 * tangent vectors etc. They ensure that 1, 2, and 3D objects are
 * always represented in 3D space, while 4 and 5D objects (needed in
 * kinetic solvers) use NDIM space.
 */
  template <unsigned NDIM> class CDIM {}; // fall-back case
  template <> class CDIM<1> { static const unsigned N = 3; };
  template <> class CDIM<2> { static const unsigned N = 3; };
  template <> class CDIM<3> { static const unsigned N = 3; };
  template <> class CDIM<4> { static const unsigned N = 4; };
  template <> class CDIM<5> { static const unsigned N = 5; };
}

#endif // LC_CDIM_H
