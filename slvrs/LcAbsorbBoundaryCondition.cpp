/**
 * @file	LcAbsorbBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAbsorbBoundaryCondition.h>

namespace Lucee
{
  static const unsigned LC_LOWER_EDGE = 0;
  static const unsigned LC_UPPER_EDGE = 1;

// set costructor name
  const char *AbsorbBoundaryCondition::id = "Absorb";

  void
  AbsorbBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::BoundaryCondition::readInput(tbl);
// ensure exactly 3 components are specified
    if (this->numComponents() != 3)
      throw Lucee::Except(
        "AbsorbBoundaryCondition::readInput: Absorbing BCs can be applied only to 3-component vector");
  }

  void
  AbsorbBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
// rotate vector to local coordinate system
    double vel[3], normVel[3];
    for (unsigned i=0; i<3; ++i)
      vel[i] = qin[this->component(i)];
    c.rotateVecToLocal(vel, normVel);

// flip sign of normal component if it going back into the domain
    const unsigned edge = getEdge();
    if ((edge == LC_LOWER_EDGE && normVel[0] > 0.)
        || (edge == LC_UPPER_EDGE && normVel[0] < 0.))
    normVel[0] = -normVel[0];

// rotate back to global frame
    c.rotateVecToGlobal(normVel, vel);

    for (unsigned i=0; i<3; ++i)
      qbc[this->component(i)] = vel[i];
  }
}
