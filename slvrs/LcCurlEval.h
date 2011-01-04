/**
 * @file	LcCurlEval.h
 *
 * @brief	Compute curl in rectangular coordinates.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */
#ifndef LC_CURL_EVAL_H
#define LC_CURL_EVAL_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/**
 * Class to evaluate curl of a vector field in DIR direction. The
 * vector field must live on a rectangular coordinate system.
 */
  template <unsigned DIR> struct CurlEval;

/** Curl in X direction */
  template <> struct CurlEval<0>
  {
/**
 * Compute curl in X-direction. Computes
 *
 * cv[n] = cv[n] + delta*curl(v) 
 *
 * where only X derivatives are used. 
 *
 * @param delta spacing factor (generally dt/dx).
 * @param v Vector field.
 * @param vr Vector feild on right of v.
 * @param cv Updated field.
 */
      static void
      eval(double delta, const double v[3], const double vr[3], double cv[3])
      {
        cv[1] = cv[1] - delta*(vr[2]-v[2]);
        cv[2] = cv[2] + delta*(vr[1]-v[1]);
      }
  };

/** Curl in Y direction */
  template <> struct CurlEval<1>
  {
/**
 * Compute curl in Y-direction. Computes
 *
 * cv[n] = cv[n] + delta*curl(v) 
 *
 * where only Y derivatives are used. 
 *
 * @param delta spacing factor (generally dt/dy).
 * @param v Vector field.
 * @param vr Vector feild on right of v.
 * @param cv Updated field.
 */
      static void
      eval(double delta, const double v[3], const double vr[3], double cv[3])
      {
        cv[0] = cv[0] + delta*(vr[2]-v[2]);
        cv[2] = cv[2] - delta*(vr[0]-v[0]);
      }
  };

/** Curl in Z direction */
  template <> struct CurlEval<2>
  {
/**
 * Compute curl in Z-direction. Computes
 *
 * cv[n] = cv[n] + delta*curl(v) 
 *
 * where only Z derivatives are used. 
 *
 * @param delta spacing factor (generally dt/dz).
 * @param v Vector field.
 * @param vr Vector feild on right of v.
 * @param cv Updated field.
 */
      static void
      eval(double delta, const double v[3], const double vr[3], double cv[3])
      {
        cv[0] = cv[0] - delta*(vr[1]-v[1]);
        cv[1] = cv[1] + delta*(vr[0]-v[0]);
      }
  };
}

#endif // LC_CURL_EVAL_H
