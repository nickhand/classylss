%{
#include "ClassParams.h"
%}

namespace std {
    %template(VectorString) vector<std::string>;
};
%rename(_update) update(const ClassParams &other);
%rename(_contains) contains(const std::string& key) const;

class ClassParams {
public:    
    
    ClassParams();
    ClassParams(const std::string& param_file);
    //ClassParams(const ClassParams &other);
    //~ClassParams();
    
    void update(const ClassParams &other);

    int add(const std::string& key, const int val); 
    int add(const std::string& key, const long val); 
    int add(const std::string& key, const long long val); 
    int add(const std::string& key, const unsigned val); 
    int add(const std::string& key, const unsigned long val); 
    int add(const std::string& key, const unsigned long long val); 
    int add(const std::string& key, const float val); 
    int add(const std::string& key, const double val);
    int add(const std::string& key, const long double val); 
    int add(const std::string& key, const bool val); 
    int add(const std::string& key, const std::string val);
    int add(const std::string& key, const char* val); 

    std::vector<std::string> keys() const;
    unsigned size() const;
    void print() const; 
    std::string pop(const std::string& key); 
    bool contains(const std::string& key) const;
    
    const std::string& value(const std::string& key) const;
    
    %pythoncode %{
    @classmethod
    def from_dict(cls, d):
        toret = cls()
    
        # default CLASS parameters
        d.setdefault('P_k_max_h/Mpc',  20.)
        d.setdefault('z_max_pk', 2.0)
    
        for k in d: toret.add(k, d[k])
        return toret
    
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
        pars = self.__class__.from_dict(kwargs)
        self._update(pars)
    
    def __contains__(self, key):
        return self._contains(key)
    
    def __getitem__(self, key):
        try:
            return self.value(key)
        except:
            raise KeyError("no such key: '%s'" %key)
    
    def __setitem__(self, key, value):
        return self.add(key, value)
    
    def __iter__(self):
        return iter(self.keys()) 
    
    def values(self): 
        return [self[k] for k in self.keys()]
    
    def items(self):
        return list(zip(self.keys(), self.values()))
    
    def __len__(self):
        return self.size()
    
    def __str__(self):
        return dict(self).__str__()

    def __repr__(self):
        return dict(self).__repr__()

    @classmethod
    def from_astropy(cls, cosmo, extra={}):
        """
        Convert an astropy cosmology to a ``ClassParams`` instance
        """
        from astropy import units, cosmology
    
        pars = cls()
        pars['h'] = cosmo.h
        pars['Omega_g'] = cosmo.Ogamma0
        if cosmo.Ob0 is not None:
            pars['Omega_b'] = cosmo.Ob0
        else:
            raise ValueError("please specify a value 'Ob0' ")
        pars['Omega_cdm'] = cosmo.Om0 - cosmo.Ob0 # should be okay for now

        # handle massive neutrinos
        if cosmo.has_massive_nu: 

            # convert to eV
            m_nu = cosmo.m_nu
            if m_nu.unit != units.eV:
                m_nu = m_nu.to(units.eV)

            # from CLASS notes:
            # one more remark: if you have respectively 1,2,3 massive neutrinos, 
            # if you stick to the default value pm equal to 0.71611, designed to give m/omega of 
            # 93.14 eV, and if you want to use N_ur to get N_eff equal to 3.046 in the early universe, 
            # then you should pass here respectively 2.0328,1.0196,0.00641
            N_ur = [2.0328, 1.0196, 0.00641]
            N_massive = (m_nu > 0.).sum()
            pars['N_ur'] = (cosmo.Neff/3.046) * N_ur[N_massive-1]

            pars['N_ncdm'] = N_massive
            pars['m_ncdm'] = ", ".join([str(k.value) for k in sorted(m_nu[m_nu > 0.], reverse=True)])
        else:
            pars['N_ur'] = cosmo.Neff
            pars['N_ncdm'] = 0
            pars['m_ncdm'] = 0.
        
        # handle dark energy
        if isinstance(cosmo, cosmology.LambdaCDM):
            pars['w0_fld'] = -1.0
            pars['wa_fld'] = 0.
        elif isinstance(cosmo, cosmology.wCDM):
            pars['w0_fld'] = cosmo.w0
            pars['wa_fld'] = 0.
            pars['Omega_Lambda'] = 0. # use Omega_fld
        elif isinstance(cosmo, cosmology.w0waCDM):
            pars['w0_fld'] = cosmo.w0
            pars['wa_fld'] = cosmo.wa
            pars['Omega_Lambda'] = 0. # use Omega_fld
        else:
            cls = cosmo.__class__.__name__
            valid = ["LambdaCDM", "wCDM", "w0waCDM"]
            msg = "dark energy equation of state not recognized for class '%s'; " %cls
            msg += "valid classes: %s" %str(valid)
            raise TypeError(msg)

        # default CLASS parameters
        extra.setdefault('P_k_max_h/Mpc',  20.)
        extra.setdefault('z_max_pk', 2.0)

        # add any extra arguments
        if len(extra):
            pars.update(extra)

        return pars
    %}
};
