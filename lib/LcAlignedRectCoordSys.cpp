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
      throw Lucee::Except("AlignedRectCoordSys::rotateSymMatrixToLocal: Not implemented");
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
      throw Lucee::Except("AlignedRectCoordSys::rotateSymMatrixToGlobal: Not implemented");
    }
  }
}
