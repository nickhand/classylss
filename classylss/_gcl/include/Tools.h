#ifndef TOOLS_H
#define TOOLS_H

#include "Common.h"
#include "parray.h"

// enum the integration method for Xi_lm
namespace IntegrationMethods { enum Type {FFTLOG=0, SIMPS, TRAPZ}; }

/* Functions to compute the quantity
    \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k) 
*/

// use FFTLog -- k, Pk will be interpolated onto a log-spaced grid
parray ComputeXiLM(int l, int m, const parray& k, const parray& pk, double dlogr, double logrc, double nc,
                    const parray& r, double smoothing, IntegrationMethods::Type method=IntegrationMethods::FFTLOG);

// use FFTlog, given log-spaced input k, Pk -- fills r, xi 
void ComputeXiLM_fftlog(int l, int m, int N, const double k[], const double pk[], 
                        double r[], double xi[], double smoothing=0.);

/*  Compute the correlation function xi(r) from a power spectrum P(k), this is 
  just Xi_0^2 in the notation above*/
parray pk_to_xi(int ell, const parray& k, const parray& pk, const parray& r, 
                double smoothing, IntegrationMethods::Type method=IntegrationMethods::FFTLOG);

/* Compute the power spectrum P(k) from a correlation function xi(r), sampled
* at logarithmically spaced points r[i]. */
parray xi_to_pk(int ell, const parray& r, const parray& xi, const parray& k, 
                double smoothing, IntegrationMethods::Type method=IntegrationMethods::FFTLOG);

#endif // TOOLS_H
