/**
 * @file	LcMathLib.cpp
 *
 * @brief	Math library for use in Lucee.
 */

// lucee includes
#include <LcMathLib.h>

#ifdef HAVE_GSL
# include <gsl/gsl_math.h>
# include <gsl/gsl_sf_legendre.h>
#endif

namespace Lucee
{
  void
  gauleg(int n, double x1, double x2, Lucee::Vector<double>& xv, Lucee::Vector<double>& wv)
  {
    const double EPS = 3.0e-15;
#ifdef HAVE_GSL
    const double PI = M_PI;
#else
    const double PI = 3.14159265358979323846;
#endif

// code below assumes arrays start from 1
    double *x = &xv[0] - 1;
    double *w = &wv[0] - 1;
// adapted from Numerical Recipies in C book
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);
    for (i = 1; i <= m; i++)
    {
      z = cos(PI*(i-0.25)/(n+0.5));
      do
      {
        p1 = 1.0;
        p2 = 0.0;
        for (j = 1; j <= n; j++)
        {
          p3 = p2;
          p2 = p1;
          p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        }
        pp = n*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1-p1/pp;
      } 
      while( fabs(z-z1) > EPS );
      x[i] = xm-xl*z;
      x[n+1-i] = xm+xl*z;
      w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
      w[n+1-i] = w[i];
    }
  }

  void
  legendre(int l, int m, const Lucee::Vector<double>& x, Lucee::Vector<double>& plm)
  {
#ifdef HAVE_GSL
    const double PI = M_PI;
    int i, j = plm.getLower(0);
    double fact = 2*sqrt(PI)/sqrt(2.0*l+1);
// use GSL routine to do actual computation
    for(i = x.getLower(0); i<x.getUpper(0); ++i)
      plm[j++] = fact*gsl_sf_legendre_sphPlm(l, m, x[i]);
#else
    throw Lucee::Except("Lucee::legendre requires linking with the GNU Scientific Library!");
#endif
  }

  double
  legendre(int l, int m, double x)
  {
#ifdef HAVE_GSL
    const double PI = M_PI;
    double fact = 2*sqrt(PI)/sqrt(2.0*l+1);
// use GSL routine to do actual computation
    return fact*gsl_sf_legendre_sphPlm(l, m, x);
#else
    throw Lucee::Except("Lucee::legendre requires linking with the GNU Scientific Library!");
#endif
  }

  double
  legendrePoly(int n, double x)
  {
#ifdef HAVE_GSL
// use GSL to do actual computation
    return gsl_sf_legendre_Pl(n, x);
#else
// should eventually provide an alternative implementation based on
// recursion relations for Legendre polynomials
    throw Lucee::Except("Lucee::legendrePoly requires linking with the GNU Scientific Library!");
#endif
  }

  double
  legendrePolyDeriv(int n, double x)
  {
    double dp0 = 0;
    double dp1 = 1;
    double dpn;

    if (n==0) return dp0;
    if (n==1) return dp1;

    dpn = ((n+1)*x*legendrePoly(n,x) - (n+1)*legendrePoly(n+1,x))/(1-x*x);

    return dpn;
  }
}
