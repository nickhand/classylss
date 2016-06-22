%{
#include "Cosmology.h"
%}

%rename(_update) update(const ClassParams &other);

%rename(transfers) TransferFit;

%nodefaultctor TransferFit;
%nodefaultdtor TransferFit;
struct TransferFit {
    enum Type {CLASS=0, EH, EH_NoWiggle, BBKS};
};


class Cosmology : public ClassEngine {

public:    
    
    Cosmology(bool verbose=false);
    Cosmology(TransferFit::Type tf, bool verbose=false);
    Cosmology(const std::string& param_file, bool verbose=false);
    Cosmology(const std::string& param_file, TransferFit::Type tf, bool verbose=false);
    Cosmology(const ClassParams& pars, bool verbose=false);
    Cosmology(const ClassParams& pars, TransferFit::Type tf, bool verbose=false);
    ~Cosmology();
    Cosmology(const Cosmology &other);

    void verbose(bool verbose=true);
    
    void SetTransferFunction(TransferFit::Type tf);
    void NormalizeTransferFunction(double sigma8);
    void SetSigma8(double sigma8);
    void update(const ClassParams& newpars);
    
    // parameter accessors
    double A_s() const;
    double ln_1e10_A_s() const;
    double delta_H() const;
    double sigma8() const;    
    TransferFit::Type GetTransferFit() const;
    const std::string& GetParamFile() const;
    parray GetDiscreteK() const;
    parray GetDiscreteTk() const;
    
    double EvaluateTransfer(double k) const;
    
      
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