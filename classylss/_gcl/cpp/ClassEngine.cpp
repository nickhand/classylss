#include "ClassEngine.h"
#include "DiscreteQuad.h"
#include "Spline.h"

#include <sstream>
#include <string>
#include <cmath>
#include <cstdio>

using namespace std;
using namespace Common;


ClassEngine::ClassEngine(bool verbose) : ready(false), verbose_(verbose) {}



void ClassEngine::compute(const ClassParams& pars) 
{    
    // initialize
    Initialize(pars);
    
    // run
    Run();
}

void ClassEngine::compute(const string & param_file)
{
    string fname(param_file);
    const ClassParams parameters(fname);

    // and initialize
    Initialize(parameters);
    
    // run
    Run();
}

void ClassEngine::compute()
{
    const ClassParams parameters;

    // and initialize
    Initialize(parameters);
    
    // run
    Run();
}

void ClassEngine::Run()
{
      // run
      int status = RunCLASS();
      if (status)
          throw_error("Error running CLASS", __FILE__, __LINE__);
      
      // store Omega0_m and Omegar_0 at z = 0 so we don't keep computing it
      Omega0_m_ = BackgroundValue(0., ba.index_bg_Omega_m);
      Omega0_r_ = BackgroundValue(0., ba.index_bg_Omega_r);
  
      ready = true;  
}

void ClassEngine::Initialize(const ClassParams& pars)
{        
    ClassParams pars_ = pars; 
    if (verbose_) {
        pars_.add("background_verbose", "1");
        pars_.add("thermodynamics_verbose", "1");
        pars_.add("perturbations_verbose", "1");
        pars_.add("transfer_verbose", "1");
        pars_.add("primordial_verbose", "1");
        pars_.add("spectra_verbose", "1");
        pars_.add("nonlinear_verbose", "1");
        pars_.add("lensing_verbose", "1");
    }
    
    // initialize the file content
    fc.size = 0;
    auto pfc_input = &fc;
            
    // prepare fp structure
    size_t n = pars_.size();
    char file_not_needed[] = "not_needed_file.dat";
    parser_init(pfc_input, n, file_not_needed, errmsg_);
  
    // set up the input 
    int i = 0;
    for (ClassParams::const_iterator iter = pars_.begin(); iter != pars_.end(); iter++) {        
        string key(iter->first);
        string val(pars_.value(iter->first));
        strcpy(pfc_input->name[i], key.c_str());
        strcpy(pfc_input->value[i], val.c_str());
        i++;
    }
    
    // initialize the input
    if (input_init(pfc_input, &pr, &ba, &th, &pt, &tr, &pm, &sp, &nl, &le, &op, errmsg_) == _FAILURE_) 
        throw invalid_argument(errmsg_);
}

// destructor
ClassEngine::~ClassEngine()
{
    Clean();
}

void ClassEngine::Clean() {
    
    if (ready) {
        Free();
        ready = false;
    }
}

// print content of file_content
void ClassEngine::PrintFC() const {
    
    printf("FILE_CONTENT SIZE=%d\n",fc.size);
    for (int i=0; i<fc.size; i++) 
        printf("%d : %s = %s\n",i,fc.name[i],fc.value[i]);
}

// main CLASS function from main.c
int ClassEngine::class_main(struct file_content *pfc,
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
			                    ErrorMsg errmsg) {
    

    if (input_init(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,pop,errmsg) == _FAILURE_) {
        info("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
        ready = false;
        return _FAILURE_;
    }

    if (background_init(ppr,pba) == _FAILURE_) {
        info("\n\nError running background_init \n=>%s\n",pba->error_message);
        ready = false;
        return _FAILURE_;
    }

    if (thermodynamics_init(ppr,pba,pth) == _FAILURE_) {
        info("\n\nError in thermodynamics_init \n=>%s\n",pth->error_message);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (perturb_init(ppr,pba,pth,ppt) == _FAILURE_) {
        info("\n\nError in perturb_init \n=>%s\n",ppt->error_message);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (primordial_init(ppr,ppt,ppm) == _FAILURE_) {
        info("\n\nError in primordial_init \n=>%s\n",ppm->error_message);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (nonlinear_init(ppr,pba,pth,ppt,ppm,pnl) == _FAILURE_)  {
        info("\n\nError in nonlinear_init \n=>%s\n",pnl->error_message);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (transfer_init(ppr,pba,pth,ppt,pnl,ptr) == _FAILURE_) {
        info("\n\nError in transfer_init \n=>%s\n",ptr->error_message);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (spectra_init(ppr,pba,ppt,ppm,pnl,ptr,psp) == _FAILURE_) {
        info("\n\nError in spectra_init \n=>%s\n",psp->error_message);
        transfer_free(&tr);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    if (lensing_init(ppr,ppt,psp,pnl,ple) == _FAILURE_) {
        info("\n\nError in lensing_init \n=>%s\n",ple->error_message);
        info("Calling spectra free #1\n");
        spectra_free(&sp);
        transfer_free(&tr);
        nonlinear_free(&nl);
        primordial_free(&pm);
        perturb_free(&pt);
        thermodynamics_free(&th);
        background_free(&ba);
        ready = false;
        return _FAILURE_;
    }

    ready = true;
    return _SUCCESS_;
}

// run CLASS to do the desired work
int ClassEngine::RunCLASS() {

    int status = this->class_main(&fc, &pr, &ba, &th, &pt, &tr, &pm, &sp, &nl, &le, &op, errmsg_);
    return status;
}

int ClassEngine::Free() {
  
    if (lensing_free(&le) == _FAILURE_) {
        throw_error(le.error_message, __FILE__, __LINE__);
    }
  
    if (nonlinear_free(&nl) == _FAILURE_) {
        throw_error(nl.error_message, __FILE__, __LINE__);
    }
  
    if (spectra_free(&sp) == _FAILURE_) {
        throw_error(sp.error_message, __FILE__, __LINE__);
    }
    
    if (primordial_free(&pm) == _FAILURE_) {
        throw_error(pm.error_message, __FILE__, __LINE__);
    }
    
    if (transfer_free(&tr) == _FAILURE_) {
        throw_error(tr.error_message, __FILE__, __LINE__);
    }

    if (perturb_free(&pt) == _FAILURE_) {
        throw_error(pt.error_message, __FILE__, __LINE__);
    }

    if (thermodynamics_free(&th) == _FAILURE_) {
        throw_error(th.error_message, __FILE__, __LINE__);
    } 

    if (background_free(&ba) == _FAILURE_) {
        throw_error(ba.error_message, __FILE__, __LINE__);
    }
    return _SUCCESS_;
}

int ClassEngine::check_le_cls_type(Cls::Type t)
{
    switch(t) {
        case Cls::TT:
            if (le.has_tt==_FALSE_) 
                throw invalid_argument("no lensed ClTT available");
            return le.index_lt_tt;
        case Cls::TE:
            if (le.has_te==_FALSE_) 
                throw invalid_argument("no lensed ClTE available"); 
            return le.index_lt_te;
        case Cls::EE:
            if (le.has_ee==_FALSE_)
                throw invalid_argument("no lensed ClEE available");
            return le.index_lt_ee;
        case Cls::BB:
            if (le.has_bb==_FALSE_)
                throw invalid_argument("no lensed ClBB available");
            return le.index_lt_bb;
        case Cls::PP:
            if (le.has_pp==_FALSE_)
                throw invalid_argument("no lensed ClPhi-Phi available");
            return le.index_lt_pp;
        case Cls::TP:
            if (le.has_tp==_TRUE_)
                throw invalid_argument("no lensed ClT-Phi available");
            return le.index_lt_tp;
    } 
}

int ClassEngine::check_sp_cls_type(Cls::Type t)
{
    switch(t) {
        case Cls::TT:
            if (sp.has_tt==_FALSE_) 
                throw invalid_argument("no ClTT available");
            return sp.index_ct_tt;
        case Cls::TE:
            if (sp.has_te==_FALSE_) 
                throw invalid_argument("no ClTE available"); 
            return sp.index_ct_te;
        case Cls::EE:
            if (sp.has_ee==_FALSE_)
                throw invalid_argument("no ClEE available");
            return sp.index_ct_ee;
        case Cls::BB:
            if (sp.has_bb==_FALSE_)
                throw invalid_argument("no ClBB available");
            return sp.index_ct_bb;
        case Cls::PP:
            if (sp.has_pp==_FALSE_)
                throw invalid_argument("no ClPhi-Phi available");
            return sp.index_ct_pp;
        case Cls::TP:
            if (sp.has_tp==_TRUE_)
                throw invalid_argument("no ClT-Phi available");
            return sp.index_ct_tp;
    } 
}

int ClassEngine::lmax() const {
    if (!ready) throw_error("`lmax` set only available after calling compute()", __FILE__, __LINE__);
    return min(sp.l_max_tot, pt.l_scalar_max);
}

int ClassEngine::lensed_lmax() const {
    if (!ready) throw_error("`lensed_lmax` set only available after calling compute()", __FILE__, __LINE__);
    return min(le.l_lensed_max, pt.l_scalar_max);
}

/*----------------------------------------------------------------------------*/
/* Functions for getting spectra */
/*----------------------------------------------------------------------------*/


parray ClassEngine::GetRawCls(const parray& ell, Cls::Type t) 
{
    if (!ready) throw_error("run compute() before accessing results", __FILE__, __LINE__);    
    int index = check_sp_cls_type(t);
    
    // check bounds
    if (ell.max() > sp.l_max_tot) {
        ostringstream msg;
        msg << "maximum ell value should be below ell = " << sp.l_max_tot;
        throw_error(msg.str(), __FILE__, __LINE__);
    }
    if (ell.min() < 2)
        throw_error("minimum ell value is ell=2", __FILE__, __LINE__);
    
    double *rcl = new double[sp.ct_size]();
    
    // quantities for tensor modes
    double **cl_md = new double*[sp.md_size];
    for (int i = 0; i < sp.md_size; ++i)
        cl_md[i] = new double[sp.ct_size]();
    
    // quantities for isocurvature modes
    double **cl_md_ic = new double*[sp.md_size];
    for (int i = 0; i < sp.md_size; ++i)
        cl_md_ic[i] = new double[sp.ct_size*sp.ic_ic_size[i]]();
    
    double tomuk = 1e6*Tcmb();
    double tomuk2 = tomuk*tomuk;
        
    parray toret = parray::zeros(ell.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < toret.size(); i++)
    {
        if (spectra_cl_at_l(&sp, ell[i], rcl, cl_md, cl_md_ic) == _FAILURE_)
            throw invalid_argument(sp.error_message);
        toret[i] = tomuk2*rcl[index];
    }
    
    delete [] rcl;
    return toret;
}

parray ClassEngine::GetLensedCls(const parray& ell, Cls::Type t)
{
    if (!ready) throw_error("run compute() before accessing results", __FILE__, __LINE__);
    int index = check_le_cls_type(t);
    
    // check bounds
    if (ell.max() > le.l_lensed_max) {
        ostringstream msg;
        msg << "maximum lensed ell value should be below ell = " << le.l_lensed_max;
        throw_error(msg.str(), __FILE__, __LINE__);
    }
    if (ell.min() < 2)
        throw_error("minimum ell value is ell=2", __FILE__, __LINE__);

    double *lcl = new double[le.lt_size]();

    double tomuk = 1e6*Tcmb();
    double tomuk2 = tomuk*tomuk;

    // return array
    parray toret = parray::zeros(ell.size());

    #pragma omp parallel for
    for (size_t i = 0; i < toret.size(); i++)
    {
        if (lensing_cl_at_l(&le, ell[i], lcl) == _FAILURE_)
            throw invalid_argument(le.error_message);
        toret[i] = tomuk2*lcl[index];
    }

    delete [] lcl;
    return toret;
}

// generic private class for computing either P_lin(k) or P_nl(k)
parray ClassEngine::GetPk(const parray& k, double z, Pktype method) {
    
    double *pk_ic = new double[sp.ic_ic_size[sp.index_md_scalars]]();
    int status;
    double thisPk;
    double h3 = pow3(h());
    parray toret(k.size());
    parray scaled_k = k*h();
    
    // check for k > kmax
    if (scaled_k.max() > k_max())
        throw_error("Computing P(k) for k > kmax", __FILE__, __LINE__);
    
    // don't use openmp --> bad things happen
    for (size_t i = 0; i < toret.size(); i++) 
    {
        // make sure to put k from h/Mpc to 1/Mpc
        if (method == Pk_linear)
            status = spectra_pk_at_k_and_z(&ba, &pm, &sp, scaled_k[i], z, &thisPk, pk_ic);
        else
            status = spectra_pk_nl_at_k_and_z(&ba, &pm, &sp, scaled_k[i], z, &thisPk);

        if (status == _FAILURE_)
            throw invalid_argument(sp.error_message);
        
        toret[i] = thisPk*h3; // put into units of (Mpc/h^3)
    }
    return toret; 
}

// compute the k, linear Pk in units of h/Mpc, (Mpc/h)^3
parray ClassEngine::GetPklin(const parray& k, double z) {
    
    if (!ready) throw_error("run compute() before accessing results", __FILE__, __LINE__);    
    if (pt.has_pk_matter == _FALSE_)
        throw_error("cannot compute linear power if matter P(k) not requested", __FILE__, __LINE__);
    return GetPk(k, z, Pk_linear);
}

// compute the k, nonlinear Pk in units of h/Mpc, (Mpc/h)^3
parray ClassEngine::GetPknl(const parray& k, double z) {
        
    if (!ready) throw_error("run compute() before accessing results", __FILE__, __LINE__);    
    if (nl.method == nl_none)
        throw_error("cannot compute nonlinear P(k) if `non linear` parameter is not set to `halofit`", __FILE__, __LINE__);
    return GetPk(k, z, Pk_nonlinear);
}

// return the k, transfer function (in CAMB format) in units of h/Mpc, unitless
parray ClassEngine::GetTk(const parray& k, double z) 
{
    if (!ready) throw_error("run compute() before accessing results", __FILE__, __LINE__);    
    if (pt.has_density_transfers != _TRUE_)
        throw_error("no density transfers were requested", __FILE__, __LINE__);
    
    int index_md = sp.index_md_scalars;
        
    // power spectrum
    parray Tk;
    if (nl.method == nl_none)
        Tk = GetPklin(k, z);
    else
        Tk = GetPknl(k, z);
    
    for (size_t i = 0; i < k.size(); i++) {
    
        // scale by the primodial power spectrum
        double pk_primordial_k;
        double thisk = k[i]*h();
        primordial_spectrum_at_k(&pm, index_md, linear, thisk, &pk_primordial_k);
        Tk[i] /= (pk_primordial_k*k[i]*pow3(h()));
        Tk[i] = sqrt(Tk[i]);
    }
    return Tk;
}

/*----------------------------------------------------------------------------*/
/* Functions for computing cosmological quantities as a function of z */
/*----------------------------------------------------------------------------*/
double ClassEngine::BackgroundValue(double z, int index) const {
    
    double tau;
    double toret;
    int index2;
    parray pvecback(ba.bg_size);
    
    // find the conformal time for this redshift
    background_tau_of_z(&ba, z, &tau);

    // call to fill pvecback
    background_at_tau(&ba, tau, ba.long_info, ba.inter_normal, &index2, pvecback);
    if (index > 0) 
        toret = pvecback[index];
    else
        toret = -1;
    return toret;
}

/*----------------------------------------------------------------------------*/

// the power spectrum normalization, sigma8
double ClassEngine::Sigma8_z(double z) const {
    

    double sigma8 = 0.;
    BackgroundValue(z, -1);
    spectra_sigma(&ba, &pm, &sp, 8./ba.h, z, &sigma8);

    return sigma8;
}

/*----------------------------------------------------------------------------*/
double ClassEngine::rho_crit(bool cgs) const {
    
    double H = 100. * (Constants::km / Constants::second / Constants::Mpc);
    double rho_crit_cgs = 3*pow(H, 2)/(8*M_PI*Constants::G);
    
    if (cgs)
        return rho_crit_cgs;
    else
        return rho_crit_cgs / (Constants::M_sun / pow(Constants::Mpc, 3));
}

/*----------------------------------------------------------------------------*/
double ClassEngine::rho_crit_z(double z, bool cgs) const {
    
    double Ez = H_z(z) / H0();
    return rho_crit(cgs) * pow(Ez, 2);
}


/*----------------------------------------------------------------------------*/
double ClassEngine::Omega_m_z(double z) const {
    
    double Ez = H_z(z) / H0();
    return Omega0_m() * pow(1.+z, 3) / pow(Ez, 2);
}

/*----------------------------------------------------------------------------*/
double ClassEngine::rho_bar_z(double z, bool cgs) const {
    
    return Omega_m_z(z) * rho_crit_z(z, cgs);
}

/*----------------------------------------------------------------------------*/
double ClassEngine::Dm_z(double z) const {
    double Dc = Dc_z(z);
    double toret(0.);
    if (ba.sgnK == 0)
        toret = Dc;
    else if (ba.sgnK == 1)
        toret = sin(sqrt(ba.K)*Dc)/sqrt(ba.K);
    else if (ba.sgnK == -1)
        toret = sinh(sqrt(-ba.K)*Dc)/sqrt(-ba.K);
    return toret;
}

parray ClassEngine::Dm_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::Dm_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::f_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::f_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::H_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::H_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::Da_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::Da_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::Dc_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::Dc_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::D_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::D_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::Sigma8_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::Sigma8_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::Omega_m_z(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::Omega_m_z;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

parray ClassEngine::rho_bar_z(const parray& z, bool cgs) const {
    double (ClassEngine::*pf)(double, bool) const = &ClassEngine::rho_bar_z;
    return EvaluateMany(bind(pf, this, placeholders::_1, cgs), z);
}

parray ClassEngine::rho_crit_z(const parray& z, bool cgs) const {
    double (ClassEngine::*pf)(double, bool) const = &ClassEngine::rho_crit_z;;
    return EvaluateMany(bind(pf, this, placeholders::_1, cgs), z);
}

double ClassEngine::dV(double z) const {
    double c_kms = Constants::c_light/(Constants::km/Constants::second);
    return c_kms * pow2((1+z)*Da_z(z)) / H_z(z) / pow3(1e3);  
}

parray ClassEngine::dV(const parray& z) const {
    double (ClassEngine::*pf)(double) const = &ClassEngine::dV;
    return EvaluateMany(bind(pf, this, placeholders::_1), z);
}

double ClassEngine::V(double zmin, double zmax, int Nz) const {
    parray z = parray::linspace(zmin, zmax, Nz);
    parray integrand = dV(z);
    return 4*M_PI*SimpsIntegrate(z, integrand);
}

parray ClassEngine::V(const parray& zmin, const parray& zmax, int Nz) const {
    
    int N = (int) zmin.size();
    parray toret(N);
    #pragma omp parallel for
    for(int i = 0; i < N; i++)
        toret[i] = V(zmin[i], zmax[i], Nz);
    return toret;
}


