/**
 * @file	lcfftwtest.cxx
 *
 * @brief	Test for FFTW usage from lucee
 */

// lucee includes
#include <LcTest.h>

// FFTW includes
# include <fftw3.h>

// std includes
#include <complex>
#include <cmath>

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcfftwtest");

  int N = 128;
  std::vector<std::complex<double> > in(N), out(N);
  std::vector<std::complex<double> > uprime_t(N), uprime(N);

// create alias to data in vector to allow passing to FFTW
  fftw_complex *f_in = reinterpret_cast<fftw_complex*> (&in[0]);
  fftw_complex *f_out = reinterpret_cast<fftw_complex*> (&out[0]);

  fftw_complex *b_in = reinterpret_cast<fftw_complex*> (&uprime_t[0]);
  fftw_complex *b_out = reinterpret_cast<fftw_complex*> (&uprime[0]);

// compute forward and backward plans
  fftw_plan fp = fftw_plan_dft_1d(N, f_in, f_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bp = fftw_plan_dft_1d(N, b_in, b_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  double L = 20.0; // domain is [-L/2, L/2]
  std::vector<double> X(N+1);
  double dx = L/N;
  for (unsigned i=0; i<N+1; ++i)
    X[i] = -L/2+i*dx;
  
// compute function (ignore final point due to periodicity)
  for (unsigned i=0; i<N; ++i)
    in[i] = 1/std::cosh(X[i]);

// compute FFT of function
  fftw_execute(fp);

// compute wave-number array with appropriate ordering
  double mpi = 3.141592654;
  std::vector<double> kx(N);
  for (unsigned i=0; i<N/2; ++i)
    kx[i] = 2*mpi/L*i;
  for (unsigned i=N/2; i>0; --i)
    kx[N-i] = -2*mpi/L*i;

// compute FFT of derivative of function
  std::complex<double> j1(0, 1);
  for (unsigned i=0; i<N; ++i)
    uprime_t[i] = j1*kx[i]*out[i];

// now compute inverse FFT
  fftw_execute(bp);

// output derivatives (division by N is needed as FFTW transforms are
// not normalized)
   for (unsigned i=0; i<N; ++i)
     std::cout << uprime[i].real()/N << std::endl;
  
  fftw_destroy_plan(fp); fftw_destroy_plan(bp);
  LC_END_TESTS;
}
