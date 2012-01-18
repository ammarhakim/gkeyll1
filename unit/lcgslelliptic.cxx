/**
 * @file	lcgslelliptic.cxx
 *
 * @brief	Test for elliptic functions as computed from GSL
 */

// lucee includes
#include <LcTest.h>

#ifdef HAVE_GSL
# include <gsl/gsl_math.h>
# include <gsl/gsl_sf_ellint.h>
#endif

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcgslelliptic");

  double x0=0.0, x1=0.99;
  unsigned n=10;
  double dx = (x1-x0)/n;

  std::cout << "Elliptic E" << std::endl;
  for (unsigned i=0; i<=n; ++i)
    std::cout << i << " " << x0+i*dx << " " 
              << gsl_sf_ellint_Ecomp(x0+i*dx, GSL_PREC_DOUBLE) << std::endl;

  std::cout << "Elliptic K" << std::endl;
  for (unsigned i=0; i<=n; ++i)
    std::cout << i << " " << x0+i*dx << " " 
              << gsl_sf_ellint_Kcomp(x0+i*dx, GSL_PREC_DOUBLE) << std::endl;

  LC_END_TESTS;
}
