%{
#include "PowerSpectrum.h"
%}

class PowerSpectrum {
public:
    PowerSpectrum();
    virtual ~PowerSpectrum();

    // translated to __call__ -> calls Evaluate(K)
    double operator()(const double k) const;
    
    // translated to __call__ -> calls EvaluateMany(K)
    parray operator()(const parray& k) const;
     
    // mass variance sigma(R)
    double Sigma(double R) const;
    parray Sigma(const parray& R) const;
    
    // 1D velocity dispersion
    double VelocityDispersion() const;
    
    // 1 / 1D velocity disp
    double NonlinearScale() const;
    
    // sigma_v as a function of k
    parray VelocityDispersion(const parray& k, double factor = 0.5) const;
    double VelocityDispersion(const double k, double factor = 0.5) const;
};

