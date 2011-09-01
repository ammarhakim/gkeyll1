/**
 * @file	lcfftwtest.cxx
 *
 * @brief	Test for FFTW usage from lucee
 */

// lucee includes
#include <LcTest.h>

#ifdef HAVE_FFTW
# include <fftw.h>
#endif

int
main(void)
{
  LC_BEGIN_TESTS("lcfftwtest");

  LC_END_TESTS;
}
