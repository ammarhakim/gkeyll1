/**
 * @file	LcMathLib.h
 *
 * @brief	Math library for use in Lucee.
 *
 * @version	$Id: LcMathLib.h 305 2010-03-01 23:08:07Z a.hakim777 $
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
 * @param x Ordinates for Gaussian quadrature.
 * @param w Weights for Gaussian quadrature.
 */
  void gauleg(int n, double x1, double x2, Lucee::Vector<double>& x, 
    Lucee::Vector<double>& w);
}

#endif // LC_MATH_LIB_H
