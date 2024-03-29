/**
 * @file	LcAlignedRectCoordSys.h
 *
 * @brief	Rectangular coordinate system.
 */

#ifndef LC_ALIGNED_RECT_COORD_SYS_H
#define LC_ALIGNED_RECT_COORD_SYS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRectCoordSys.h>

namespace Lucee
{
/**
 * A rectangular coordinate system obtained from a 90-degree rotation.
 */
  class AlignedRectCoordSys : public Lucee::RectCoordSys
  {
    public:
/**
 * Create a coordinate system such that the local X-axis point along
 * 'dir' direction of the global coordinate system. The other two axis
 * are aligned to create a right-handed coordinate system. Note that
 * there is still an ambiguity in the other direction.
 *
 * @param dir Local X-axis points along 'dir' direction.
 */
      AlignedRectCoordSys(unsigned dir);

/**
 * Rotate vector to local coordinate system (defined by this object)
 * from global coordinate system.
 *
 * @param inVec Vector to rotate.
 * @param outVec Rotated vector.
 */
      void rotateVecToLocal(const double inVec[3], double outVec[3]) const;

/**
 * Rotate vector from local coordinate system (defined by this object)
 * to global coordinate system.
 *
 * @param inVec Vector to rotate.
 * @param outVec Rotated vector.
 */
      void rotateVecToGlobal(const double inVec[3], double outVec[3]) const;

/**
 * Rotate symmetric 3x3 matrix to local coordinate system (defined by
 * this object) from global coordinate system. The matrix is stored as
 * an array of 6 components [A11, A12, A13, A22, A23, A33].
 *
 * @param inSM Symmetric matrix to rotate.
 * @param outSM Rotated matrix.
 */
      void rotateSymMatrixToLocal(const double inSM[6], double outSM[6]) const;

/**
 * Rotate symmetric 3x3 matrix from local coordinate system (defined
 * by this object) to global coordinate system.  The matrix is stored
 * as an array of 6 components [A11, A12, A13, A22, A23, A33].
 *
 * @param inSM Symmetric matrix to rotate.
 * @param outSM Rotated matrix.
 */
      void rotateSymMatrixToGlobal(const double inSM[6], double outSM[6]) const;

/**
 * Rotate symmetric 3x3x3 matrix to local coordinate system (defined by
 * this object) from global coordinate system. The matrix is stored as
 * an array of 10 components [A111, A112, A113, A122, A123, A133, A222, A223, A233, A333].
 *
 * @param inSM Symmetric tensor to rotate.
 * @param outSM Rotated tensor.
 */
      void rotateSymTensorToLocal(const double inSM[10], double outSM[10]) const;

/**
 * Rotate symmetric 3x3x3 tensor from local coordinate system (defined
 * by this object) to global coordinate system.  The matrix is stored
 * as an array of 10 components [A111, A112, A113, A122, A123, A133, A222, A223, A233, A333].
 *
 * @param inSM Symmetric tensor to rotate.
 * @param outSM Rotated tensor.
 */
      void rotateSymTensorToGlobal(const double inSM[10], double outSM[10]) const;


/**
 * Returns the alignment direction dir for this object
 */
      int getAlignmentDirection() const;

    private:
/** Alignment direction */
      unsigned dir;
  };
}

#endif // LC_ALIGNED_RECT_COORD_SYS_H
