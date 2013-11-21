/**
 * @file	LcMathPhysConstants.h
 *
 * @brief	Mathematical and physical constants.
 */

#ifndef LC_MATH_PHYS_CONSTANTS_H
#define LC_MATH_PHYS_CONSTANTS_H

namespace Lucee
{
// math constants
  static const double PI = 3.14159265358979323846264338328;
  static const double E = 2.71828182845904523536028747135;
  static const double EULER = 0.57721566490153286060651209008;

// physical constants
  static const double SPEED_OF_LIGHT = 2.99792458e8; // m / s
  static const double PLANCKS_CONSTANT_H = 6.62606896e-34; // kg m^2 / s
  static const double ELECTRON_MASS = 9.10938188e-31; // kg
  static const double PROTON_MASS = 1.672621777e-27; // kg
  static const double ELEMENTARY_CHARGE = 1.602176565e-19; // C
  static const double BOLTZMANN_CONSTANT = 1.3806488e-23; // J/K
  static const double EPSILON0 = 8.854187817e-12; // F/m
  static const double MU0 = 12.566370614e-7; // H/m
  static const double EV2KELVIN = 1.1604519e4; // eV/K
}

#endif // LC_MATH_PHYS_CONSTANTS_H
