%{
#include "LinearPS.h"
%}

%rename(_update) update(const ClassParams &other);

class LinearPS : public PowerSpectrum {
public:
    LinearPS(const Cosmology& cosmo, double z=0);
    
    const double& GetSigma8AtZ() const;
    const Cosmology& GetCosmology() const;
    const double& GetRedshift() const;
    
    void SetSigma8AtZ(double sigma8_z);
    
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
    
    void update(const ClassParams& newpars);
    
    %pythoncode %{
            
        def update(self, *args, **kwargs):
            if len(args):
                if len(args) != 1:
                    raise ValueError("only one positional argument, a dictionary")
                d = args[0]
                if not isinstance(d, dict):
                    raise TypeError("first argument must be a dictionary")
                kwargs.update(d)
                
            if not len(kwargs):
                raise ValueError("no parameters provided to update")
            pars = ClassParams.from_dict(kwargs)
            self._update(pars)
    %}
};

