%{
#include "ClassEngine.h"
%}

namespace std {
    %template(VectorUnsigned) std::vector<unsigned>;
};

%rename(_compute) compute(const std::string& param_file);
%rename(_compute) compute();
%rename(_compute) compute(const ClassParams& pars);
%rename(_clean) Clean();

%pythoncode %{
import contextlib
%}

%rename(cltypes) Cls;

%nodefaultctor Cls;
%nodefaultdtor Cls;
struct Cls {
    enum Type {TT=0,EE,TE,BB,PP,TP};
};

class ClassEngine  {

public:

    static std::string Alpha_inf_hyrec_file;
    static std::string R_inf_hyrec_file;
    static std::string two_photon_tables_hyrec_file;
    static std::string sBBN_file;
    
    ClassEngine(bool verbose=false);  
    ~ClassEngine();

    // compute from a set of parameters
    void compute(const std::string& param_file);
    void compute();
    void compute(const ClassParams& pars);
    
    // is ready?
    bool isready() const;
    
    // set verbosity
    void verbose(bool verbose=true);
    
    // Cls functions
    parray GetRawCls(const parray& ell, Cls::Type cl=Cls::TT);
    parray GetLensedCls(const parray& ell, Cls::Type cl=Cls::TT);
    
    // matter density/transfer
    parray GetPklin(const parray& k, double z);
    parray GetPknl(const parray& k, double z);
    parray GetTk(const parray& k, double z);
    
    int lmax() const;
    int lensed_lmax() const;
    
    double H0() const;
    double h() const;
    double Tcmb() const;
    double Omega0_b() const;
    double Omega0_cdm() const;
    double Omega0_ur() const;
    double Omega0_m() const;
    double Omega0_r() const;
    double Omega0_g() const;
    double Omega0_lambda() const;
    double Omega0_fld() const;
    double Omega0_k() const;
    double w0_fld() const;
    double wa_fld() const;
    double n_s() const;
    double k_pivot() const;
    double A_s() const;
    double ln_1e10_A_s() const;
    double sigma8() const;
    double k_max() const;
    double k_min() const;
    double z_drag() const;
    double rs_drag() const;
    double tau_reio() const;
    double z_reio() const;
    
    double f_z(double z) const;
    parray f_z(const parray& z) const;
    
    double H_z(double z) const;
    parray H_z(const parray& z) const;
    
    double Da_z(double z) const;
    parray Da_z(const parray& z) const;

    double Dc_z(double z) const;
    parray Dc_z(const parray& z) const;
    
    double Dm_z(double z) const;
    parray Dm_z(const parray& z) const;
    
    double D_z(double z) const;
    parray D_z(const parray& z) const;
    
    double Sigma8_z(double z) const;
    parray Sigma8_z(const parray& z) const;
    
    double Omega_m_z(double z) const;
    parray Omega_m_z(const parray& z) const;
    
    double rho_bar_z(double z, bool cgs = false) const;
    parray rho_bar_z(const parray& z, bool cgs = false) const;
    
    double rho_crit_z(double z, bool cgs = false) const;
    parray rho_crit_z(const parray& z, bool cgs = false) const;
    
    double dV(double z) const;
    parray dV(const parray& z) const;
    
    double V(double zmin, double zmax, int Nz=1024) const;
    parray V(const parray& zmin, const parray& zmax, int Nz=1024) const;

    void Clean();
    
    %pythoncode %{

        @contextlib.contextmanager
        def compute(self, *args):
            
            try:
                self._compute(*args)
                yield
            except Exception as e:
                raise e
            finally:
                self._clean()
    %}
};
