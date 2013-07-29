/**
 * @file	LcZeroTangentBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcZeroTangentBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *ZeroTangentBoundaryCondition::id = "ZeroTangent";

  void
  ZeroTangentBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::BoundaryCondition::readInput(tbl);
// ensure exactly 3 components are specified
    if (this->numComponents() != 3)
      throw Lucee::Except(
        "ZeroTangentBoundaryCondition::readInput: Zero-tangent BCs can be applied only to 3-component vector");
  }

  void
  ZeroTangentBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// rotate vector to local coordinate system
    double vel[3], normVel[3];
    for (unsigned i=0; i<3; ++i)
      vel[i] = qin[this->component(i)];
    c.rotateVecToLocal(vel, normVel);

// flip sign of tangent components
    normVel[1] = -normVel[1];
    normVel[2] = -normVel[2];

// rotate back to global frame
    c.rotateVecToGlobal(normVel, vel);

    for (unsigned i=0; i<3; ++i)
      qbc[this->component(i)] = vel[i];
  }
}
