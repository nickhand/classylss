/*------------------------------------------------------------------------------
Description:
    class ClassEngine : encapsulation of class calls


Author List:
    Stephane Plaszczynski (plaszczy@lal.in2p3.fr)
    Nick Hand (nhand@berkeley.edu)

History (add to end):
    creation:   ven. nov. 4 11:02:20 CET 2011 
    modified: nov 21 14
-----------------------------------------------------------------------------*/
#ifndef CLASSENGINE_H
#define CLASSENGINE_H


#include "ClassParams.h"
#include "Common.h"
#include "parray.h"
#include <functional>

namespace Cls { 
  enum Type {TT=0,EE,TE,BB,PP,TP}; //P stands for phi (lensing potential)
}

/*----------------------------------------------------------------------------*/
/* class to serve as the engine for CLASS */
/*----------------------------------------------------------------------------*/
class ClassEngine 
{
public:
    
    ClassEngine(bool verbose=false);  
    ~ClassEngine();
    
    // compute from a set of parameters
    void compute(const std::string& param_file);
    void compute();
    void compute(const ClassParams& pars);
    
    // is ready? i.e, whether compute has been called
    bool isready() const { return ready; }
    
    // set verbosity
    inline void verbose(bool verbose=true) { verbose_=verbose; }

    // return the primary Cls
    parray GetRawCls(const parray& ell, Cls::Type t=Cls::TT);
    
    // return the lensing Cls
    parray GetLensedCls(const parray& ell, Cls::Type t=Cls::TT);
    
    // linear/nonlinear matter power spectra
    parray GetPklin(const parray& k, double z);
    parray GetPknl(const parray& k, double z);

    // get the total density transfer function
    parray GetTk(const parray& k, double z);

    // parameter accessors
    int lmax() const; 
    int lensed_lmax() const;

    /*------------------------------------------------------------------------*/
    /* present day quantities */
    /*------------------------------------------------------------------------*/
    // the present-day Hubble constant in units of km/s/Mpc
    inline double H0() const { return ba.H0*Constants::c_light/(Constants::km/Constants::second); }
    // the dimensionless Hubble constant
    inline double h() const { return ba.h; }
    // CMB temperature today in Kelvin
    inline double Tcmb() const { return ba.T_cmb; }
    // present baryon density parameter (\Omega_b)
    inline double Omega0_b() const { return ba.Omega0_b; }
    // present cold dark matter density fraction (\Omega_cdm)
    inline double Omega0_cdm() const { return ba.Omega0_cdm; }
    // present ultra-relativistic neutrino density fraction (\Omega_{\nu, ur})
    inline double Omega0_ur() const { return ba.Omega0_ur; }
    // present non-relativistic density fraction (\Omega_b + \Omega_cdm + \Omega_{\nu nr}) 
    inline double Omega0_m() const { return Omega0_m_; }
    // relativistic density fraction (\Omega_{\gamma} + \Omega_{\nu r})
    inline double Omega0_r() const { return Omega0_r_; }
    // photon density parameter today
    inline double Omega0_g() const { return ba.Omega0_g; }
    // cosmological constant density fraction
    inline double Omega0_lambda() const { return ba.Omega0_lambda; }
    // dark energy fluid density fraction (valid if Omega0_lambda is unspecified)
    inline double Omega0_fld() const { return ba.Omega0_fld; }
    // curvature density fraction
    inline double Omega0_k() const { return ba.Omega0_k; }
    // current fluid equation of state parameter
    inline double w0_fld() const { return ba.w0_fld; }
    // fluid equation of state parameter derivative
    inline double wa_fld() const { return ba.wa_fld; }
    
    // the spectral index of the primordial power spectrum
    inline double n_s() const { return pm.n_s; }
    // pivot scale in Mpc-1 
    inline double k_pivot() const { return pm.k_pivot; }
    // scalar amplitude = curvature power spectrum at pivot scale
    inline double A_s() const { return pm.A_s; }
    // convenience function returns log (1e10*A_s)
    inline double ln_1e10_A_s() const { return log(1e10*A_s()); }
    // convenience function to return sigma8 at z = 0
    inline double sigma8() const { return sp.sigma8; }
    // maximum k value computed in h/Mpc
    inline double k_max() const { return 0.995*exp(sp.ln_k[sp.ln_k_size-1])/h(); }
    // minim k value computed in h/Mpc
    inline double k_min() const { return 1.005*exp(sp.ln_k[0])/h(); }
    
    // baryon drag redshift
    inline double z_drag() const { return th.z_d; }
    // comoving sound horizon at baryon drag
    inline double rs_drag() const { return th.rs_d; } 
    // reionization optical depth
    inline double tau_reio() const { return th.tau_reio; }
    // reionization redshift
    inline double z_reio() const { return th.z_reio; }
    
    // critical density at z = 0 in units of h^2 M_sun / Mpc^3 if cgs = False, 
    // or in units of h^2 g / cm^3
    double rho_crit(bool cgs = false) const; 
    
    /*------------------------------------------------------------------------*/
    /* background quantities as a function of z */
    /*------------------------------------------------------------------------*/    
    // the growth rate f(z) (unitless)
    inline double f_z(double z) const { return BackgroundValue(z, ba.index_bg_f); }
    parray f_z(const parray& z) const;
    
    // Hubble constant H(z) (km/s/Mpc)
    double H_z(double z) const { return (z != 0.) ? BackgroundValue(z, ba.index_bg_H)*Constants::c_light/(Constants::km/Constants::second) : H0(); }
    parray H_z(const parray& z) const;
        
    // angular diameter distance in Mpc -- this is Dm/(1+z)
    double Da_z(double z) const { return BackgroundValue(z, ba.index_bg_ang_distance); }
    parray Da_z(const parray& z) const;

    // conformal distance in flat case in Mpc
    double Dc_z(double z) const { return BackgroundValue(z, ba.index_bg_conf_distance); }
    parray Dc_z(const parray& z) const;
    
    // comoving radius coordinate in Mpc -- equal to conformal distance in flat case
    double Dm_z(double z) const;
    parray Dm_z(const parray& z) const;
    
    // growth function D(z) / D(0) (normalized to unity at z = 0)
    double D_z(double z) const { return (z != 0.) ? BackgroundValue(z, ba.index_bg_D) / BackgroundValue(0., ba.index_bg_D) : 1.; }
    parray D_z(const parray& z) const;
    
    // sigma8 (z) as derived from the scalar amplitude
    // may not be equal to the desired sigma8
    double Sigma8_z(double z) const;
    parray Sigma8_z(const parray& z) const;

    // Omega0_m as a function of z
    double Omega_m_z(double z) const;
    parray Omega_m_z(const parray& z) const;
    
    // mean matter density in units of h^2 M_sun / Mpc^3 if cgs = False, or
    // in units of g / cm^3
    double rho_bar_z(double z, bool cgs = false) const;
    parray rho_bar_z(const parray& z, bool cgs = false) const;
    
    // critical matter density in units of h^2 M_sun / Mpc^3 if cgs = False, or
    // in units of g / cm^3
    double rho_crit_z(double z, bool cgs = false) const;
    parray rho_crit_z(const parray& z, bool cgs = false) const;
    
    // the comoving volume element per unit solid angle per unit redshift in Gpc^3
    double dV(double z) const;
    parray dV(const parray& z) const;
    
    // the comoving volume between two redshifts (full sky)
    double V(double zmin, double zmax, int Nz=1024) const;
    parray V(const parray& zmin, const parray& zmax, int Nz=1024) const;
            
    void Clean();
    
protected:
  
    // the main class structures
    // make them mutable such that they can be changed for different
    // redshifts, etc
    mutable struct file_content fc;
    mutable struct precision pr;        /* for precision parameters */
    mutable struct background ba;       /* for cosmological background */
    mutable struct thermo th;           /* for thermodynamics */
    mutable struct perturbs pt;         /* for source functions */
    mutable struct transfers tr;        /* for transfer functions */
    mutable struct primordial pm;       /* for primordial spectra */
    mutable struct spectra sp;          /* for output spectra */
    mutable struct nonlinear nl;        /* for non-linear spectra */
    mutable struct lensing le;          /* for lensed spectra */
    mutable struct output op;           /* for output files */

    ErrorMsg errmsg_;                   /* for error messages */

    // helpers
    bool ready;
    bool verbose_;
    
    // print content of file_content (mostly for debugging)
    void PrintFC() const;
    
    
    // check whether we computed given Cls
    int check_sp_cls_type(Cls::Type t); 
    int check_le_cls_type(Cls::Type t); 

    // call once per compute() call
    void Initialize(const ClassParams& pars);
    void Run();
    int RunCLASS();
    
    // Free
    int Free();

    int class_main(struct file_content *pfc,
    	           struct precision * ppr,
    	           struct background * pba,
    	           struct thermo * pth,
    	           struct perturbs * ppt,
    	           struct transfers * ptr,
    	           struct primordial * ppm,
    	           struct spectra * psp,
    	           struct nonlinear * pnl,
    	           struct lensing * ple,
    	           struct output * pop,
    	           ErrorMsg errmsg);
        
    // nonlinear vs linear P(k)
    enum Pktype {Pk_linear=0, Pk_nonlinear};
    parray GetPk(const parray& k, double z, Pktype method = Pk_linear);
    double GetCl(Cls::Type t, const long &l); 
        
    // functions for returning cosmological quantities
    double BackgroundValue(double z, int index) const;
        
    // omega0_m and omega0_r are constants
    double Omega0_m_, Omega0_r_;
    
    // evaluate background quantities at several redshifts
    template<typename Function>
    parray EvaluateMany(Function f, const parray& z) const {
      int Nz = (int) z.size();
      parray toret(Nz);
      #pragma omp parallel for
      for(int i = 0; i < Nz; i++)
          toret[i] = f(z[i]);
      return toret;
    }
    
};

#endif

