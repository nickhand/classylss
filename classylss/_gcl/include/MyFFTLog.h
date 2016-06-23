#ifndef MYFFTLOG_H
#define MYFFTLOG_H

#include <complex>
typedef std::complex<double> dcomplex;

/* Compute the discrete Hankel transform of the function a(r).  See the FFTLog
 * documentation (or the Fortran routine of the same name in the FFTLog
 * sources) for a description of exactly what this function computes.
 * If u is NULL, the transform coefficients will be computed anew and discarded
 * afterwards.  If you plan on performing many consecutive transforms, it is
 * more efficient to pre-compute the u coefficients. */
void fht(int N, const double r[], const dcomplex a[], double k[], dcomplex b[], double mu,
         double q = 0, double kcrc = 1, bool noring = true, dcomplex* u = NULL);

/* Pre-compute the coefficients that appear in the FFTLog implementation of
 * the discrete Hankel transform.  The parameters N, mu, and q here are the
 * same as for the function fht().  The parameter L is defined (for whatever
 * reason) to be N times the logarithmic spacing of the input array, i.e.
 *   L = N * log(r[N-1]/r[0])/(N-1) */
void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, dcomplex u[]);


#endif // MYFFTLOG_H