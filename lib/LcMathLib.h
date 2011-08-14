/**
 * @file	LcMathLib.h
 *
 * @brief	Math library for use in Lucee.
 */

#ifndef LC_MATH_LIB_H
#define LC_MATH_LIB_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcVector.h>
#include <LcVec3.h>

// std includes
#include <algorithm>
#include <cmath>

namespace Lucee
{
/**
 * Compute ordinates and weights for n-th order Gaussian quadrature in
 * the interval [x1, x2].
 *
 * @param n Number of ordinates.
 * @param x1 Left end of interval.
 * @param x2 Right end of interval.
 * @param x On output, ordinates for Gaussian quadrature.
 * @param w On output, weights for Gaussian quadrature.
 */
  void gauleg(int n, double x1, double x2, 
    Lucee::Vector<double>& x, Lucee::Vector<double>& w);

/**
 * Compute normalized assosiated Legendre functions. These are defined by
 *
 * P_l^m = \sqrt{(l-m)!/(l+m)!} (1-x)^{m/2} d^
 *
 * @param l Coefficient in P_l^m
 * @param m Coefficient in P_l^m
 * @param x Values at which P_l^m are required.
 * @param plm On output, P_l^m(x[i]).
 */
  void legendre(int l, int m, const Lucee::Vector<double>& x, Lucee::Vector<double>& plm);

/**
 * Compute normalized assosiated Legendre function. These are defined by
 *
 * P_l^m = \sqrt{(l-m)!/(l+m)!} (1-x)^{m/2} d^
 *
 * @param l Coefficient in P_l^m
 * @param m Coefficient in P_l^m
 * @param x Value at which P_l^m are required.
 * @return plm P_l^m(x).
 */
  double legendre(int l, int m, double x);

/**
 * Compute minimum of three numbers.
 *
 * @param a First of three numbers.
 * @param b Second of three numbers.
 * @param c Third of three numbers.
 * @return minimum of three numbers.
 */
  template <typename T>
  const T& min3(const T& a, const T& b, const T& c)
  {
    return std::min<T>(a, std::min<T>(b, c));
  }

/**
 * Compute maximum of three numbers.
 *
 * @param a First of three numbers.
 * @param b Second of three numbers.
 * @param c Third of three numbers.
 * @return maximum of three numbers.
 */
  template <typename T>
  const T& max3(const T& a, const T& b, const T& c)
  {
    return std::max<T>(a, std::max<T>(b, c));
  }

/**
 * Compute area of a triangle.
 *
 * @param a Vertex of triangle.
 * @param b Vertex of triangle.
 * @param c Vertex of triangle.
 * @return area of triangle.
 */
  template <typename T>
  T calcTriArea(const Lucee::Vec3<T>& a, const Lucee::Vec3<T>& b, const Lucee::Vec3<T>& c)
  {
    Lucee::Vec3<T> s1 = b-a;
    Lucee::Vec3<T> s2 = c-a;
    return 0.5*s1.cross(s2).getNorm();
  }

/**
 * Compute area of a quadrilateral (a,b,c,d), with sides (a,b), (b,c),
 * (c,d) and (d,a).
 *
 * @param a Vertex of quadrilateral.
 * @param b Vertex of quadrilateral.
 * @param c Vertex of quadrilateral.
 * @param d Vertex of quadrilateral.
 * @return area of quadrilateral.
 */
  template <typename T>
  T calcQuadArea(const Lucee::Vec3<T>& a, const Lucee::Vec3<T>& b, 
    const Lucee::Vec3<T>& c, const Lucee::Vec3<T>& d)
  {
    Lucee::Vec3<T> s1 = b-a;
    Lucee::Vec3<T> d1 = c-a;
    Lucee::Vec3<T> s2 = d-a;
    return 0.5*(s1.cross(d1).getNorm() + d1.cross(s2).getNorm());
  }

/**
 * Check if number is nan.
 *
 * @param x number to check.
 * @return true if is nan, false otherwise.
 */
  template <typename FLT>
  bool isNan(const FLT& x)
  { // This might need to be rewritten for Windows as it does not support C99 standard
    return std::isnan(x);
  }
}

#endif // LC_MATH_LIB_H
