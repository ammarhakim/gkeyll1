/**
 * @file	LcZeroNormalBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcZeroNormalBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *ZeroNormalBoundaryCondition::id = "ZeroNormal";

  void
  ZeroNormalBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::BoundaryCondition::readInput(tbl);
// ensure exactly 3 components are specified
    if (this->numComponents() != 3)
      throw Lucee::Except(
        "ZeroNormalBoundaryCondition::readInput: Zero-normal BCs can be applied only to 3-component vector");
    useExtrapolation = false;
    if (tbl.hasBool("useExtrapolation"))
      useExtrapolation = tbl.getBool("useExtrapolation");
  }

  void
  ZeroNormalBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& qin1, const Lucee::ConstFieldPtr<double>& qin,
    Lucee::FieldPtr<double>& qbc)
  {
// rotate vector to local coordinate system
    double vel[3], normVel[3];
    if (useExtrapolation)
      for (unsigned i=0; i<3; ++i)
        vel[i] = 1.5*qin[this->component(i)]-0.5*qin1[this->component(i)];
    else
      for (unsigned i=0; i<3; ++i)
        vel[i] = qin[this->component(i)];
    c.rotateVecToLocal(vel, normVel);

// flip sign of normal component
    normVel[0] = -normVel[0];

// rotate back to global frame
    c.rotateVecToGlobal(normVel, vel);

    for (unsigned i=0; i<3; ++i)
      qbc[this->component(i)] = vel[i];
  }
}
