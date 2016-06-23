#include "ZeldovichPS.h"
#include "LinearPS.h"
#include "FortranFFTLog.h"

using namespace Common;

static const int NUM_PTS = 1024;
static const double RMIN = 1e-2;
static const double RMAX = 1e5;
static const int NMAX = 15;


/*----------------------------------------------------------------------------*/
/* Zeldovich base class */
/*----------------------------------------------------------------------------*/
ZeldovichPS::ZeldovichPS(const Cosmology& C_, double z_, bool approx_lowk_) 
    : C(C_), sigma8_z(C_.Sigma8_z(z_)), Plin(C_, z_), approx_lowk(approx_lowk_), k0_low(5e-3)
{    
    // initialize the R array
    InitializeR();

    // compute integrals at z = 0 first
    LinearPS P_L(C_, 0.);
    double norm = pow2(sigma8_z/C.sigma8());
    
    // the integrals we need
    sigma_sq = norm*P_L.VelocityDispersion();    
    X0 = norm*P_L.X_Zel(r);
    XX = X0 + 2.*sigma_sq;
    YY = norm*P_L.Y_Zel(r); 
}


ZeldovichPS::ZeldovichPS(const Cosmology& C_, bool approx_lowk_, double sigma8_z_, double k0_low_,
                         double sigmasq, const parray& X0_, const parray& XX_, const parray& YY_) 
                        : C(C_), sigma8_z(sigma8_z_), Plin(C_, 0.), approx_lowk(approx_lowk_), 
                         k0_low(k0_low_), sigma_sq(sigmasq), X0(X0_), XX(XX_), YY(YY_)
{
    InitializeR();
    Plin.SetSigma8AtZ(sigma8_z);      
}

void ZeldovichPS::InitializeR() 
{
    
    r = parray(NUM_PTS);
    nc = 0.5*double(NUM_PTS+1);
    double logrmin = log10(RMIN); 
    double logrmax = log10(RMAX);
    logrc = 0.5*(logrmin+logrmax);
    dlogr = (logrmax - logrmin)/NUM_PTS; 
    
    for (int i = 1; i <= NUM_PTS; i++)
        r[i-1] = pow(10., (logrc+(i-nc)*dlogr));
}


ZeldovichPS::~ZeldovichPS() {}

static void nearest_interp_1d(int nd, double *xd, double *yd, int ni, double *xi, double *yi) 
{
    double d, d2;
    int i, j, k, l;
    
    for (i = 0; i < ni; i++) {
        k = 0;
        d = fabs(xi[i] - xd[k]);
        for(j = 1; j < nd; j++) {
            d2 = fabs(xi[i] - xd[j]);
            if (d2 < d) {
                k = j;
                d = d2;
            }
        }
        if (xi[i] > xd[k]) {
            l = k+1;
            yi[i] = yd[l] + (yd[k] - yd[l])/(xd[k] - xd[l])*(xi[i]-xd[l]);
        } else if (xi[i] < xd[k]) {
            l = k-1;
            yi[i] = yd[l] + (yd[k] - yd[l] )/(xd[k] - xd[l])*(xi[i]-xd[l]);
        } else
            yi[i] = yd[k];
    }
}

parray ZeldovichPS::EvaluateMany(const parray& k) const 
{    
    int n = (int)k.size();
    parray pk(n);
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
        pk[i] = Evaluate(k[i]);
    }
    return pk;
}

void ZeldovichPS::SetSigma8AtZ(double new_sigma8_z) 
{     
    double ratio = pow2(new_sigma8_z / sigma8_z);
    // integrals are proportional to the square of the ratio of sigma8
    X0 *= ratio;
    XX *= ratio;
    YY *= ratio;
    sigma_sq *= ratio;
    
    // store the sigma8_z
    sigma8_z = new_sigma8_z;
    
    // set the Plin
    Plin.SetSigma8AtZ(sigma8_z);
}


double ZeldovichPS::fftlog_compute(double k, double factor) const 
{    
    double q = 0; // unbiased
    double mu;
    
    parray a(NUM_PTS);
    parray kmesh(NUM_PTS);

    // logspaced between RMIN and RMAX
    double this_Pk = 0.;
    double toadd;   
    for (int n = 0; n <= NMAX; n++) {
        
        // the order of the Bessel function
        mu = 0.5 + double(n);
        
        // compute a(r)
        if (n == 0) 
            Fprim(a, r, k);
        else
            Fsec(a, r, k, double(n));
        
        // do the fft
        FortranFFTLog fftlogger(NUM_PTS, dlogr*log(10.), mu, q, 1.0, 1);
        
        bool ok = fftlogger.Transform(a, 1);
        if (!ok) error("FFTLog failed\n");
        double kr = fftlogger.KR();
        double logkc = log10(kr) - logrc;
        
        for (int j = 1; j <= NUM_PTS; j++) 
            kmesh[j-1] = pow(10., (logkc+(j-nc)*dlogr));
             
        // sum it up
        double out;
        nearest_interp_1d(NUM_PTS, (double*)(kmesh), (double*)(a), 1, &k, &out);
        toadd = factor*sqrt(0.5*M_PI)*pow(k, -1.5)*out; 
        
        this_Pk += toadd;   
        if (fabs(toadd/this_Pk) < 0.005) break;
    }
        
    return this_Pk;
}

double ZeldovichPS::Evaluate(double k) const 
{
    if (k >= k0_low || !approx_lowk)
        return fftlog_compute(k, 4*M_PI);
    else
        return LowKApprox(k);
}

void ZeldovichPS::Fprim(parray& a, const parray& r, double k) const 
{    
    for (int i = 0; i < NUM_PTS; i++) {
        a[i] = pow(r[i], 1.5) * (exp(-0.5*pow2(k)*(XX[i] + YY[i])) - exp(-pow2(k)*sigma_sq));
    }
    
}

void ZeldovichPS::Fsec(parray& a, const parray& r, double k, double n) const 
{    
    for (int i = 0; i < NUM_PTS; i++) {
        a[i] = pow(r[i], 1.5-n)*pow(k*YY[i], n)*exp(-0.5*pow2(k)*(XX[i] + YY[i]));
    }
}

double ZeldovichPS::LowKApprox(double k) const 
{    
    return (1 - pow2(k)*sigma_sq + 0.5*pow4(k)*pow2(sigma_sq))*Plin(k) + 0.5*Plin.Q3_Zel(k);
}






 
