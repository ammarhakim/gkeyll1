/**
 * @file	LcIdealInnerPlanetBoundaryCondition.cpp
 *
 * @brief	Class for applying BCs to inner planet surface based on ideal-MHD
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcIdealInnerPlanetBoundaryCondition.h>

namespace Lucee
{
// set costructor name
  const char *IdealInnerPlanetBoundaryCondition::id = "IdealInnerPlanet";

  void
  IdealInnerPlanetBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::BoundaryCondition::readInput(tbl);
// tell base class all components will be specified    
    std::vector<unsigned> cd(18);
    for (unsigned i=0; i<18; ++i)
      cd[i] = i;
    this->setComponents(cd);
  }

  void
  IdealInnerPlanetBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c,
    const Lucee::ConstFieldPtr<double>& qin1, const Lucee::ConstFieldPtr<double>& qin,
    Lucee::FieldPtr<double>& qbc)
  {
// rotate vector to local coordinate system
    double vel[3], normVel[3], velGhost[3];
    const Lucee::ConstFieldPtr<double> *staticEB = this->getExtraInputField();

// first float everything
    for (unsigned i=0; i<18; ++i)
      qbc[i] = qin[i];

// now compute average bulk velocity
    for (unsigned i=0; i<3; ++i)
      vel[i] = (qin[1+i]+qin[6+i])/(qin[0+i]+qin[5+i]);
    c.rotateVecToLocal(vel, normVel);
// flip sign of normal component
    normVel[0] = -normVel[0];
// rotate back to global frame
    c.rotateVecToGlobal(normVel, velGhost);

// compute total magnetic field
    double B[3];
    for (unsigned i=0; i<3; ++i)
      B[i] = (*staticEB)[3+i]+qin[13+i];

// now compute ghost electric field
    qbc[10] = qin[10] + (vel[1]-velGhost[1])*B[2] - (vel[2]-velGhost[2])*B[1];
    qbc[11] = qin[11] + (vel[2]-velGhost[2])*B[0] - (vel[0]-velGhost[0])*B[2];
    qbc[12] = qin[12] + (vel[0]-velGhost[0])*B[1] - (vel[1]-velGhost[1])*B[0];
  }
}
