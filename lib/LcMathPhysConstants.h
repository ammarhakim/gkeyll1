/**
 * @file	LcMathPhysConstants.h
 *
 * @brief	Mathematical and physical constants.
 */

#ifndef LC_MATH_PHYS_CONSTANTS_H
#define LC_MATH_PHYS_CONSTANTS_H

// boost includes
#include <boost/math/constants/constants.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/codata_constants.hpp>

// make extracting value from Boost constants easier
#define _LCGBV(v) v.value().value()

namespace Lucee
{
// math constants
  static const double PI = boost::math::double_constants::pi;
  static const double E = boost::math::double_constants::e;
  static const double EULER = boost::math::double_constants::euler;

// physical constants
  static const double SPEED_OF_LIGHT = _LCGBV(boost::units::si::constants::codata::c);
  static const double PLANCKS_CONSTANT_H = _LCGBV(boost::units::si::constants::codata::h);
  static const double ELECTRON_MASS = _LCGBV(boost::units::si::constants::codata::m_e);
  static const double PROTON_MASS = _LCGBV(boost::units::si::constants::codata::m_p);
  static const double ELEMENTARY_CHARGE = _LCGBV(boost::units::si::constants::codata::e);
  static const double BOLTZMANN_CONSTANT = 1.3806488e-23;
  static const double EPSILON0 = _LCGBV(boost::units::si::constants::codata::epsilon_0);
  static const double MU0 = _LCGBV(boost::units::si::constants::codata::mu_0);
  static const double EV2KELVIN = ELEMENTARY_CHARGE/BOLTZMANN_CONSTANT;
}

#endif // LC_MATH_PHYS_CONSTANTS_H
