%{
#include "LinearPS.h"
%}

%rename(_update) update(const ClassParams &other);

class LinearPS : public PowerSpectrum {
public:
    LinearPS(const Cosmology& cosmo, double z=0);
    ~LinearPS();
    
    double Evaluate(double k) const;
    
    const double& GetSigma8AtZ() const;
    const Cosmology& GetCosmology() const;
    const double& GetRedshift() const;
    
    void SetSigma8AtZ(double sigma8_z);
        
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

