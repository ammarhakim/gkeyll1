/**
 * @file	LcRectCoordSys.cpp
 *
 * @brief	Rectangular corrdinate system.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcRectCoordSys.h>

namespace Lucee
{
  void
  RectCoordSys::fillWithUnitVecs(double xu[3], double yu[3], double zu[3]) const
  {
    for (unsigned i=0; i<3; ++i)
    {
      xu[i] = xunit[i];
      yu[i] = yunit[i];
      zu[i] = zunit[i];
    }
  }

  void
  RectCoordSys::fillWithUnitVecs(Lucee::Vec3& xu, Lucee::Vec3& yu, Lucee::Vec3& zu) const
  {
    for (unsigned i=0; i<3; ++i)
    {
      xu[i] = xunit[i];
      yu[i] = yunit[i];
      zu[i] = zunit[i];
    }
  }

  RectCoordSys::RectCoordSys()
  {
    for (unsigned i=0; i<3; ++i)
    {
      xunit[i] = 0.0;
      yunit[i] = 0.0;
      zunit[i] = 0.0;
    }
    xunit[0] = yunit[1] = zunit[2] = 1.0;
  }

  RectCoordSys::RectCoordSys(double xu[3], double yu[3], double zu[3])
  {
    if (checkVecs(xu, yu, zu) == false)
      throw Lucee::Except
        ("RectCoordSys::RectCoordSys: Suppled vectors must be orthonormal");
    for (unsigned i=0; i<3; ++i)
    {
      xunit[i] = xu[i];
      yunit[i] = yu[i];
      zunit[i] = zu[i];
    }
  }

  void
  RectCoordSys::setUnitVecs(double xu[3], double yu[3], double zu[3])
  {
    if (checkVecs(xu, yu, zu) == false)
      throw Lucee::Except
        ("RectCoordSys::setUnitVecs: Suppled vectors must be orthonormal");
    for (unsigned i=0; i<3; ++i)
    {
      xunit[i] = xu[i];
      yunit[i] = yu[i];
      zunit[i] = zu[i];
    }    
  }

  bool
  RectCoordSys::checkVecs(double xu[3], double yu[3], double zu[3])
  {
    return true;
  }
}
