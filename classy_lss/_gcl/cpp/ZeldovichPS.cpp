#include "ZeldovichPS.h"
#include "LinearPS.h"

using namespace Common;

static const int NUM_PTS = 1024;
static const double RMIN = 1e-2;
static const double RMAX = 1e5;
static const int NMAX = 15;


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


void ZeldovichPS::InitializeR() 
{
    r = parray::logspace(RMIN, RMAX, NUM_PTS);
}


ZeldovichPS::~ZeldovichPS() {}

parray ZeldovichPS::EvaluateMany(const parray& k) const 
{    
    int n = (int)k.size();
    parray pk(n);
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
    
    dcomplex* a = new dcomplex[NUM_PTS];
    dcomplex* b = new dcomplex[NUM_PTS];
    double* kmesh = new double[NUM_PTS];
    double* b_real = new double[NUM_PTS];

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
        
        fht(NUM_PTS, (const double*)(r), a, kmesh, b, mu, q, 1., true, NULL);
        
        // spline it
        for (int j = 0; j < NUM_PTS; j++)
            b_real[j] = b[j].real();
        Spline spl = LinearSpline(NUM_PTS, kmesh, b_real);
     
        toadd = factor*sqrt(0.5*M_PI)*pow(k, -1.5)*spl(k);        
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

void ZeldovichPS::Fprim(dcomplex a[], const double r[], double k) const
{    
    for (int i = 0; i < NUM_PTS; i++) {
        a[i] = pow(r[i], 1.5) * (exp(-0.5*pow2(k)*(XX[i] + YY[i])) - exp(-pow2(k)*sigma_sq));
    }
    
}

void ZeldovichPS::Fsec(dcomplex a[], const double r[], double k, double n) const
{    
    for (int i = 0; i < NUM_PTS; i++) {
        a[i] = pow(r[i], 1.5-n)*pow(k*YY[i], n)*exp(-0.5*pow2(k)*(XX[i] + YY[i]));
    }
}

double ZeldovichPS::LowKApprox(double k) const 
{    
    return (1 - pow2(k)*sigma_sq + 0.5*pow4(k)*pow2(sigma_sq))*Plin(k) + 0.5*Plin.Q3_Zel(k);
}






 
