/**
 * @file	LcRectCoordSys.h
 *
 * @brief	Rectangular coordinate system.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_RECT_COORD_SYS_H
#define LC_RECT_COORD_SYS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcVec3.h>

namespace Lucee
{
/**
 * Represents a rectangular coordinate system by three mutually
 * orthogonal unit vectors.
 */
  class RectCoordSys
  {
    public:
/**
 * Destroy object.
 */
      virtual ~RectCoordSys();

/**
 * Get coordinate system defined by this object. The returned vectors
 * are of unit length and form a (xu, yu, zu) forms a right-handed
 * coordinate system.
 *
 * @param xu On output, unit vector in X-direction.
 * @param yu On output, unit vector in Y-direction.
 * @param zu On output, unit vector in Z-direction.
 */
      void fillWithUnitVecs(double xu[3], double yu[3], double zu[3]) const;

/**
 * Get coordinate system defined by this object. The returned vectors
 * are of unit length and form a (xu, yu, zu) forms a right-handed
 * coordinate system.
 *
 * @param xu On output, unit vector in X-direction.
 * @param yu On output, unit vector in Y-direction.
 * @param zu On output, unit vector in Z-direction.
 */
      void fillWithUnitVecs(Lucee::Vec3<double>& xu, Lucee::Vec3<double>& yu, Lucee::Vec3<double>& zu) const;

/**
 * Rotate vector to local coordinate system (defined by this object)
 * from global coordinate system.
 *
 * @param inVec Vector to rotate.
 * @param outVec Rotated vector.
 */
      virtual void rotateVecToLocal(const double inVec[3], double outVec[3]) const = 0;

/**
 * Rotate vector from local coordinate system (defined by this object)
 * to global coordinate system.
 *
 * @param inVec Vector to rotate.
 * @param outVec Rotated vector.
 */
      virtual void rotateVecToGlobal(const double inVec[3], double outVec[3]) const = 0;

    protected:
/**
 * Default ctor create coordinate system that is same as global
 * coordinate system.
 */
      RectCoordSys();

/**
 * Create a new object with specified unit vectors.
 *
 * @param xu Unit vector in X-direction.
 * @param yu Unit vector in Y-direction.
 * @param zu Unit vector in Z-direction.
 */
      RectCoordSys(double xu[3], double yu[3], double zu[3]);

/**
 * Set coordinate system defined by this object. The supplied vectors
 * must be of unit length and (xu, yu, zu) must form a right-handed
 * coordinate system.
 *
 * @param xu Unit vector in X-direction.
 * @param yu Unit vector in Y-direction.
 * @param zu Unit vector in Z-direction.
 */
      void setUnitVecs(double xu[3], double yu[3], double zu[3]);

    private:
/** Unit vector in local x direction */
      double xunit[3];
/** Unit vector in local y direction */
      double yunit[3];
/** Unit vector in local z direction */
      double zunit[3];

/**
 * Check if supplied vectors are unit vectors and form a right handed
 * coordinate system.
 *
 * @param xu Vector in X-direction.
 * @param yu Vector in Y-direction.
 * @param zu Vector in Z-direction.
 * @return true if (xu,yu,zu) for orthonormal coordinate system, false otherwise.
 */
      bool checkVecs(double xu[3], double yu[3], double zu[3]);
  };
}

#endif // LC_RECT_COORD_SYS_H
