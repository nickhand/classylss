%{
#include "Cosmology.h"
%}

%rename(_update) update(const ClassParams &other);

class Cosmology : public ClassEngine {

public:    
    enum class TransferFit {CLASS, EH, EH_NoWiggle, BBKS};
    
    Cosmology(bool verbose=false);
    Cosmology(TransferFit tf, bool verbose=false);
    Cosmology(const std::string& param_file, bool verbose=false);
    Cosmology(const std::string& param_file, TransferFit tf, bool verbose=false);
    Cosmology(const ClassParams& pars, bool verbose=false);
    Cosmology(const ClassParams& pars, TransferFit tf, bool verbose=false);
    ~Cosmology();
    Cosmology(const Cosmology &other);

    void verbose(bool verbose=true);
    
    void SetTransferFunction(TransferFit tf);
    void NormalizeTransferFunction(double sigma8);
    void SetSigma8(double sigma8);
    void update(const ClassParams& newpars);
    
    // parameter accessors
    double A_s() const;
    double ln_1e10_A_s() const;
    double delta_H() const;
    double sigma8() const;    
    TransferFit GetTransferFit() const;
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

%pythoncode 
%{
    from argparse import Namespace
    transfers = Namespace()
    for k in Cosmology.__dict__:
        if k.startswith("TransferFit_"):
            setattr(transfers, k.split('TransferFit_')[1], getattr(Cosmology, k))
%}

