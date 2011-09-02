/**
 * @file	lcfftwtest.cxx
 *
 * @brief	Test for FFTW usage from lucee
 */

// lucee includes
#include <LcTest.h>

// include FFTW library stuff
#ifdef HAVE_FFTW3
# include <fftw3.h>
#endif

// std includes
#include <complex>

int
main(void)
{
  LC_BEGIN_TESTS("lcfftwtest");

  int N = 128;
  std::vector<std::complex<double> > in(N), out(N);
  fftw_plan p;

// create alias to data in vector to allow passing to FFTW
  fftw_complex *f_in = reinterpret_cast<fftw_complex*> (&in[0]);
  fftw_complex *f_out = reinterpret_cast<fftw_complex*> (&out[0]);

// compute plan
  p = fftw_plan_dft_1d(N, f_in, f_out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_destroy_plan(p);

  LC_END_TESTS;
}
