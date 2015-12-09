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
  static unsigned Q111 = 0;
  static unsigned Q112 = 1;
  static unsigned Q113 = 2;
  static unsigned Q122 = 3;
  static unsigned Q123 = 4; 
  static unsigned Q133 = 5;
  static unsigned Q222 = 6;
  static unsigned Q223 = 7;
  static unsigned Q233 = 8;
  static unsigned Q333 = 9;

  
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
      // Assume user knows that this coordinate system will only
      // be used to store a direction and nothing else.
      xu[0] = 1.0;
      yu[1] = 1.0;
      zu[2] = 1.0;
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
  
  void 
  AlignedRectCoordSys::rotateSymTensorToLocal(const double inSM [10], double outSM [10]) const
  {
    if (dir==0)
    { // no rotation
      for (unsigned i=0; i<10; ++i)
        outSM[i] = inSM[i];
    }
    else if (dir == 1)
    { 
        // outsm[Qijk] = R[im]R[jn]R[kl] insm[Qijk]
      outSM[Q111] = inSM[6];
      outSM[Q112] = -inSM[3];
      outSM[Q113] = inSM[7];
      outSM[Q122] = inSM[1];
      outSM[Q123] = -inSM[4];
      outSM[Q133] = inSM[8];
      outSM[Q222] = -inSM[0];
      outSM[Q223] = inSM[2];
      outSM[Q233] = -inSM[5];
      outSM[Q333] = inSM[9];
    } else {
      outSM[Q111] = inSM[9];
      outSM[Q112] = inSM[8];
      outSM[Q113] = -inSM[5];
      outSM[Q122] = inSM[7];
      outSM[Q123] = -inSM[4];
      outSM[Q133] = inSM[2];
      outSM[Q222] = inSM[6];
      outSM[Q223] = -inSM[3];
      outSM[Q233] = inSM[1];
      outSM[Q333] = -inSM[0];
    }
  }

  void 
  AlignedRectCoordSys::rotateSymTensorToGlobal(const double inSM [10], double outSM [10]) const
  {
    if (dir==0)
    { // no rotation
      for (unsigned i=0; i<10; ++i)
        outSM[i] = inSM[i];
    }
    else if (dir==1)
    { 
      outSM[Q111] = -inSM[6];
      outSM[Q112] = inSM[3];
      outSM[Q113] = inSM[7];
      outSM[Q122] = -inSM[1];
      outSM[Q123] = -inSM[4];
      outSM[Q133] = -inSM[8];
      outSM[Q222] = inSM[0];
      outSM[Q223] = inSM[2];
      outSM[Q233] = inSM[5];
      outSM[Q333] = inSM[9];
    } else {
      outSM[Q111] = -inSM[9];
      outSM[Q112] = inSM[8];
      outSM[Q113] = inSM[5];
      outSM[Q122] = -inSM[7];
      outSM[Q123] = -inSM[4];
      outSM[Q133] = -inSM[2];
      outSM[Q222] = inSM[6];
      outSM[Q223] = inSM[3];
      outSM[Q233] = inSM[1];
      outSM[Q333] = inSM[0];

    }
  }
  int
  AlignedRectCoordSys::getAlignmentDirection() const
  {
    return dir;
  }
}
