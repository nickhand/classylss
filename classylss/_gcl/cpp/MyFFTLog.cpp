#include <cmath>

#include "Common.h"
#include "MyFFTLog.h"
#include <fftw3.h>


/* Computes the Gamma function using the Lanczos approximation */
static dcomplex gamma(dcomplex z) {
    /* Lanczos coefficients for g = 7 */
    static double p[] = {
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
       -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
       -176.61502916214059906584551354,
        12.507343278686904814458936853,
       -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
    };

    if(z.real() < 0.5)
        return M_PI / (std::sin(M_PI*z)*gamma(1. - z));

    z -= 1;
    dcomplex x = p[0];
    for(int n = 1; n < 9; n++)
        x += p[n] / (z + double(n));
    dcomplex t = z + 7.5;
    return sqrt(2*M_PI) * std::pow(t, z+0.5) * std::exp(-t) * x;
}

static dcomplex lngamma(dcomplex z) {
    return std::log(gamma(z));
}

static void lngamma(double x, double y, double* lnr, double* arg) {
    dcomplex w = lngamma(dcomplex(x,y));
    if(lnr) *lnr = w.real();
    if(arg) *arg = w.imag();
}

static double goodkr(int N, double mu, double q, double L, double kr) {
    double xp = (mu+1+q)/2;
    double xm = (mu+1-q)/2;
    double y = M_PI*N/(2*L);
    double lnr, argm, argp;
    lngamma(xp, y, &lnr, &argp);
    lngamma(xm, y, &lnr, &argm);
    double arg = log(2/kr) * N/L + (argp + argm)/M_PI;
    double iarg = round(arg);
    if(arg != iarg)
        kr *= exp((arg - iarg)*L/N);
    return kr;
}

void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, dcomplex u[]) {
    double y = M_PI/L;
    double k0r0 = kcrc * exp(-L);
    double t = -2*y*log(k0r0/2);

    if(q == 0) {
        double x = (mu+1)/2;
        double lnr, phi;
        for(int m = 0; m <= N/2; m++) {
            lngamma(x, m*y, &lnr, &phi);
            u[m] = std::polar(1.0, m*t + 2*phi);
        }
    }
    else {
        double xp = (mu+1+q)/2;
        double xm = (mu+1-q)/2;
        double lnrp, phip, lnrm, phim;
        for(int m = 0; m <= N/2; m++) {
            lngamma(xp, m*y, &lnrp, &phip);
            lngamma(xm, m*y, &lnrm, &phim);
            u[m] = std::polar(exp(q*log(2) + lnrp - lnrm), m*t + phip - phim);
        }
    }

    for(int m = N/2+1; m < N; m++)
        u[m] = std::conj(u[N-m]);
    if((N % 2) == 0)
        u[N/2] = dcomplex(std::real(u[N/2]), 0.0);
}

void fht(int N, const double r[], const dcomplex a[], double k[], dcomplex b[], double mu,
         double q, double kcrc, bool noring, dcomplex* u)
{
    double L = log(r[N-1]/r[0]) * N/(N-1.);
    dcomplex* ulocal = NULL;
    if(u == NULL) {
        if(noring)
            kcrc = goodkr(N, mu, q, L, kcrc);
        ulocal = new dcomplex[N];
        compute_u_coefficients(N, mu, q, L, kcrc, ulocal);
        u = ulocal;
    }
    
    /* Compute the convolution b = a*u using FFTs */
    fftw_plan forward_plan = fftw_plan_dft_1d(N, (fftw_complex*) a, (fftw_complex*) b,  -1, FFTW_ESTIMATE);
    fftw_plan reverse_plan = fftw_plan_dft_1d(N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE);
    
    fftw_execute(forward_plan);
    
    for(int m = 0; m < N; m++)
        b[m] *= u[m] / double(N);       // divide by N since FFTW doesn't normalize the inverse FFT
    fftw_execute(reverse_plan);
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(reverse_plan);

    /* Reverse b array */
    dcomplex tmp;
    for(int n = 0; n < N/2; n++) {
        tmp = b[n];
        b[n] = b[N-n-1];
        b[N-n-1] = tmp;
    }

    /* Compute k's corresponding to input r's */
    double dlnr = log(r[N-1]/r[0]) / N;
    double nc = 0.5*double(N+1);
    double logkc = log(kcrc/sqrt(r[N-1]*r[0]));
    for(int j = 1; j <= N; j++) 
        k[j-1] = exp(logkc + (j-nc)*dlnr);

    delete[] ulocal;
}

void fftlog_ComputeXiLM(int l, int m, int N, const double k[], const double pk[], double r[], double xi[]) {
    dcomplex* a = new dcomplex[N];
    dcomplex* b = new dcomplex[N];

    for(int i = 0; i < N; i++)
        a[i] = pow(k[i], m - 0.5) * pk[i];
    fht(N, k, a, r, b, l + 0.5);
    for(int i = 0; i < N; i++)
        xi[i] = std::real(pow(2*M_PI*r[i], -1.5) * b[i]);

    delete[] b;
    delete[] a;
}

void pk2xi(int N, const double k[], const double pk[], double r[], double xi[]) {
    fftlog_ComputeXiLM(0, 2, N, k, pk, r, xi);
}

void xi2pk(int N, const double r[], const double xi[], double k[], double pk[]) {
    static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    fftlog_ComputeXiLM(0, 2, N, r, xi, k, pk);
    for(int j = 0; j < N; j++)
        pk[j] *= TwoPiCubed;
}