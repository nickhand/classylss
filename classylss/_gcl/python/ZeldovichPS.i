%{
#include "ZeldovichPS.h"
%}

class ZeldovichPS {
public:
    
    ZeldovichPS(const Cosmology& C, double z, bool approx_lowk=false);
    ~ZeldovichPS();
    
    // translated to __call__ -> calls Evaluate(K)
    double operator()(const double k) const;
    
    // translated to __call__ -> calls EvaluateMany(K)
    parray operator()(const parray& k) const;
    
    // functions for doing low-k approximation
    void SetLowKApprox(bool approx_lowk_=true);
    void SetLowKTransition(double k0);
    double LowKApprox(double k) const;
    
    const double& GetSigma8AtZ() const;
    const Cosmology& GetCosmology() const;
    const bool& GetApproxLowKFlag() const;
    const double& GetK0Low() const;
    parray GetXZel() const;
    parray GetYZel() const;
    parray GetX0Zel() const;
    const double& GetSigmaSq() const;
 
    void SetSigma8AtZ(double sigma8);
    
};
