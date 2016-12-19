#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <complex>

#include "Tools.h"
#include "PowerSpectrum.h"
#include "DiscreteQuad.h"
#include "Spline.h"
#include "FortranFFTLog.h"

using std::bind;
using std::cref;
using namespace std::placeholders;
using namespace Common;

void ComputeXiLM_fftlog(int l, int m, const parray& k, const parray& pk,
                            double r[], double xi[], double q, double smoothing)
{
    int N(k.size());
    parray a(N);
  
    // spacing quantities
    double nc = 0.5*double(N+1);
    double logkmin = log10(k[0]);
    double logkmax = log10(k[-1]);
    double dlogr = (logkmax-logkmin)/(N-1);
    double logrc = 0.5*(logkmin+logkmax);
    
    FortranFFTLog fftlogger(N, dlogr*log(10.), l+0.5, q, 1.0, 1);
    
    // the fftlog magic
    for(int i = 0; i < N; i++)
        a[i] = pow(k[i], m - 0.5) * pk[i] * exp(-pow2(k[i]*smoothing));
    
    bool ok = fftlogger.Transform(a, 1);
    if (!ok) 
        error("FFTLog failed in ComputeXiLM_fftlog (`ok` flag False)\n");
    
    double kr = fftlogger.KR();
    double logkc = log10(kr) - logrc;
    
    for (int j = 1; j <= N; j++) 
        r[j-1] = pow(10., (logkc+(j-nc)*dlogr));
    
    for(int i = 0; i < N; i++)
        xi[i] = std::real(pow(2*M_PI*r[i], -1.5) * a[i]);
}


parray ComputeXiLM(int l, int m, const parray& k_, const parray& pk_, const parray& r, 
                            double smoothing, IntegrationMethods::Type method) 
{
    
    // fftlog
    if (method == IntegrationMethods::FFTLOG) {
        int N(k_.size());
        double r_[N];
        double xi_[N];

        // force log-spacing using a spline
        parray k(N);
        double nc = 0.5*double(N+1);
        double logkmin = log10(k_.min()); 
        double logkmax = log10(k_.max());
        double logrc = 0.5*(logkmin+logkmax);
        double dlogr = (logkmax - logkmin)/(N-1); 
    
        for (int i = 1; i <= N; i++)
            k[i-1] = pow(10., (logrc+(i-nc)*dlogr));
        
        auto pk_spline = CubicSpline(k_, pk_);
        auto pk = pk_spline(k);

        // call fftlog on the double[] arrays
        ComputeXiLM_fftlog(l, m, k, pk, r_, xi_, 0., smoothing);

        // return the result at desired domain values using a spline
        auto spline = CubicSpline(N, r_, xi_);
        return spline(r);

    // discrete integral
    } else {
        
        // array sizes
        int Nk(pk_.size()), Nr(r.size());
        parray xi(Nr);

        // choose appropriate spherical Bessel function
        double (*sj)(double x);
        if (l == 0)      sj = SphericalBesselJ0;
        else if (l == 1) sj = SphericalBesselJ1;
        else if (l == 2) sj = SphericalBesselJ2;
        else if (l == 3) sj = SphericalBesselJ3;
        else if (l == 4) sj = SphericalBesselJ4;
        else if (l == 6) sj = SphericalBesselJ6;
        else if (l == 8) sj = SphericalBesselJ8;
        else {
            error("ComputeXiLM: l = %d not supported\n", l);
            return xi;
        }

        // integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
        #pragma omp parallel for
        for(int i = 0; i < Nr; i++) {
  
            // the integrand for this r
            parray integrand(Nk);
            for(int j = 0; j < Nk; j++) 
                integrand[j] = exp(-pow2(k_[j]*smoothing)) * pk_[j] * pow(k_[j], m) / (2*M_PI*M_PI) * sj(k_[j]*r[i]);
  
            // integrate using simpson's rule
            if (method == IntegrationMethods::SIMPS)
                xi[i] = SimpsIntegrate(k_, integrand);
            else if (method == IntegrationMethods::TRAPZ)
                xi[i] = TrapzIntegrate(k_, integrand);
            else
                error("the integration method in `ComputeXiLM` must be one of [`fftlog=0`, `simps=1`, `trapz=2`]\n");
        }
        return xi;
    }
}

parray pk_to_xi(int ell, const parray& k, const parray& pk, const parray& r, 
                    double smoothing, IntegrationMethods::Type method) 
{
    return pow(-1, ell/2)*ComputeXiLM(ell, 2, k, pk, r, smoothing, method);
}

parray xi_to_pk(int ell, const parray& r, const parray& xi, const parray& k, 
                    double smoothing, IntegrationMethods::Type method)
{
    static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    return pow(-1, ell/2) * TwoPiCubed * ComputeXiLM(ell, 2, r, xi, k, smoothing, method);
}

