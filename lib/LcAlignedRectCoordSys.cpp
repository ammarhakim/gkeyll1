/**
 * @file	LcAlignedRectCoordSys.cpp
 *
 * @brief	Rectangular coordinate system.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAlignedRectCoordSys.h>
#include <LcExcept.h>
#include <LcFixedVector.h>

namespace Lucee
{
// indexing into a symmetric 3x3 matrix (the assumption here is that
// the components of the symmetric matrix are stored in a vector of
// size 6 with entries in the order 11, 12, 13, 22, 23, 33).
  static unsigned S11 = 0;
  static unsigned S12 = 1;
  static unsigned S13 = 2;
  static unsigned S21 = 1;
  static unsigned S22 = 3;
  static unsigned S23 = 4;
  static unsigned S31 = 2;
  static unsigned S32 = 4;
  static unsigned S33 = 5;

  AlignedRectCoordSys::AlignedRectCoordSys(unsigned dir)
    : RectCoordSys(), dir(dir)
  {
    Lucee::FixedVector<3, double> xu(0.0), yu(0.0), zu(0.0);
    if (dir == 0)
    {
      xu[0] = 1.0;
      yu[1] = 1.0;
      zu[2] = 1.0;
    }
    else if (dir == 1)
    {
      xu[1] = 1.0;
      yu[0] = -1.0;
      zu[2] = 1.0;
    }
    else if (dir == 2)
    {
      xu[2] = 1.0;
      yu[1] = 1.0;
      zu[0] = -1.0;
    }
    else
    {
      throw Lucee::Except
        ("AlignedRectCoordSys::AlignedRectCoordSys: Alignment direction must be 0, 1, or 2");
    }
// call base class method to set vectors
    setUnitVecs(&xu[0], &yu[0], &zu[0]);
  }

  void
  AlignedRectCoordSys::rotateVecToLocal(const double inVec[3], double outVec[3]) const
  {
    if (dir == 0)
    {
      outVec[0] = inVec[0];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
    }
    else if (dir == 1)
    {
      outVec[0] = inVec[1];
      outVec[1] = -inVec[0];
      outVec[2] = inVec[2];
    }
    else
    {
      outVec[0] = inVec[2];
      outVec[1] = inVec[1];
      outVec[2] = -inVec[0];
    }
  }

  void
  AlignedRectCoordSys::rotateVecToGlobal(const double inVec[3], double outVec[3]) const
  {
    if (dir == 0)
    {
      outVec[0] = inVec[0];
      outVec[1] = inVec[1];
      outVec[2] = inVec[2];
    }
    else if (dir == 1)
    {
      outVec[0] = -inVec[1];
      outVec[1] = inVec[0];
      outVec[2] = inVec[2];
    }
    else
    {
      outVec[0] = -inVec[2];
      outVec[1] = inVec[1];
      outVec[2] = inVec[0];
    }
  }

  void
  AlignedRectCoordSys::rotateSymMatrixToLocal(const double inSM[6], double outSM[6]) const
  {
    if (dir==0)
    {
      for (unsigned i=0; i<6; ++i)
        outSM[i] = inSM[i];
    }
    else
    {
// first rotate all the columns
      double V0[3], V1[3], V2[3];
      double R0[3], R1[3], R2[3];
      
      V0[0] = inSM[S11];
      V0[1] = inSM[S21];
      V0[2] = inSM[S31];

      V1[0] = inSM[S12];
      V1[1] = inSM[S22];
      V1[2] = inSM[S32];

      V2[0] = inSM[S13];
      V2[1] = inSM[S23];
      V2[2] = inSM[S33];

      rotateVecToLocal(V0, R0);
      rotateVecToLocal(V1, R1);
      rotateVecToLocal(V2, R2);

// now rotate the resulting rows
      V0[0] = R0[0];
      V0[1] = R1[0];
      V0[2] = R2[0];

      V1[0] = R0[1];
      V1[1] = R1[1];
      V1[2] = R2[1];

      V2[0] = R0[2];
      V2[1] = R1[2];
      V2[2] = R2[2];

      rotateVecToLocal(V0, R0);
      rotateVecToLocal(V1, R1);
      rotateVecToLocal(V2, R2);

// copy over rotated data
      outSM[0] = R0[0];
      outSM[1] = R0[1];
      outSM[2] = R0[2];
      outSM[3] = R1[1];
      outSM[4] = R1[2];
      outSM[5] = R2[2];
    }
  }

  void
  AlignedRectCoordSys::rotateSymMatrixToGlobal(const double inSM[6], double outSM[6]) const
  {
    if (dir==0)
    {
      for (unsigned i=0; i<6; ++i)
        outSM[i] = inSM[i];
    }
    else
    {
// first rotate all the columns
      double V0[3], V1[3], V2[3];
      double R0[3], R1[3], R2[3];
      
      V0[0] = inSM[S11];
      V0[1] = inSM[S21];
      V0[2] = inSM[S31];

      V1[0] = inSM[S12];
      V1[1] = inSM[S22];
      V1[2] = inSM[S32];

      V2[0] = inSM[S13];
      V2[1] = inSM[S23];
      V2[2] = inSM[S33];

      rotateVecToGlobal(V0, R0);
      rotateVecToGlobal(V1, R1);
      rotateVecToGlobal(V2, R2);

// now rotate the resulting rows
      V0[0] = R0[0];
      V0[1] = R1[0];
      V0[2] = R2[0];

      V1[0] = R0[1];
      V1[1] = R1[1];
      V1[2] = R2[1];

      V2[0] = R0[2];
      V2[1] = R1[2];
      V2[2] = R2[2];

      rotateVecToGlobal(V0, R0);
      rotateVecToGlobal(V1, R1);
      rotateVecToGlobal(V2, R2);

// copy over rotated data
      outSM[0] = R0[0];
      outSM[1] = R0[1];
      outSM[2] = R0[2];
      outSM[3] = R1[1];
      outSM[4] = R1[2];
      outSM[5] = R2[2];
    }
  }
}
