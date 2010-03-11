/**
 * @file	LcMathLib.h
 *
 * @brief	Math library for use in Lucee.
 *
 * @version	$Id: LcMathLib.h 312 2010-03-02 18:37:24Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_MATH_LIB_H
#define LC_MATH_LIB_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcVector.h>

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
}

#endif // LC_MATH_LIB_H
