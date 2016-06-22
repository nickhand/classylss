#ifndef LINEAR_PS_H
#define LINEAR_PS_H


#include "Common.h"
#include "PowerSpectrum.h"
#include "Cosmology.h"

/// Linear power spectrum interpolated from a transfer function
class LinearPS : public PowerSpectrum {
public:
    
    LinearPS(const Cosmology& cosmo, double z=0);
    ~LinearPS();

    // evaluate for single k
    double Evaluate(double k) const;
    
    // accessors
    const double& GetRedshift() const { return z; }
    const double& GetSigma8AtZ() const { return sigma8_z; }
    const Cosmology& GetCosmology() const { return C; }
    
    // convenience to rescale
    void SetSigma8AtZ(double sigma8_z_) { sigma8_z=sigma8_z_; }
    
    // update the cosmology
    void update(const ClassParams& newpars) { C.update(newpars); }

private:
    
    Cosmology C;
    double z;
    double sigma8_z; 
};

#endif // LINEAR_PS_H
