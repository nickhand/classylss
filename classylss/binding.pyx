#cython: embedsignature=True
cimport cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset, strncpy, strdup

from classylss import get_data_files

_DATA_FILES = get_data_files()

from cclassy cimport *

DEF _Mpc_over_m_ = 3.085677581282e22  #/**< conversion factor from meters to megaparsecs */
                          #/* remark: CAMB uses 3.085678e22: good to know if you want to compare  with high accuracy */
DEF _Gyr_over_Mpc_ = 3.06601394e2 #/**< conversion factor from megaparsecs to gigayears
				  #       (c=1 units, Julian years of 365.25 days) */
DEF _c_ = 2.99792458e8          #  /**< c in m/s */
DEF _G_ = 6.67428e-11          #   /**< Newton constant in m^3/Kg/s^2 */
DEF _eV_ = 1.602176487e-19     #   /**< 1 eV expressed in J */

#/* parameters entering in Stefan-Boltzmann constant sigma_B */
DEF _k_B_ = 1.3806504e-23
DEF _h_P_ = 6.62606896e-34

DEF _MAXTITLESTRINGLENGTH_ = 8000

class ClassRuntimeError(RuntimeError):
    def __init__(self, message=""):
        self.message = message

    def __str__(self):
        return 'Class Error in Class: ' + self.message

class ClassBadValueError(ValueError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass


def val2str(val):
    if isinstance(val, (list, tuple)):
        return ','.join([str(i) for i in val])
    return str(val)

cdef int _build_file_content(pars, file_content * fc) except -1:
    fc.size = 0

    fc.filename = <char*>malloc(sizeof(FileArg))

    strncpy(fc.filename, "NOFILE", sizeof(FileArg))

    _pars = {
        "Alpha_inf hyrec file": _DATA_FILES['Alpha_inf_hyrec_file'],
        "R_inf hyrec file" : _DATA_FILES['R_inf_hyrec_file'],
        "two_photon_tables hyrec file" : _DATA_FILES['two_photon_tables_hyrec_file'],
        "sBBN file": _DATA_FILES['sBBN_file'],
        }

    _pars.update(pars)
    pars = _pars

    fc.size = len(pars)
    fc.name = <FileArg*> malloc(sizeof(FileArg)*len(pars))
    assert(fc.name!=NULL)

    fc.value = <FileArg*> malloc(sizeof(FileArg)*len(pars))
    assert(fc.value!=NULL)

    fc.read = <short*> malloc(sizeof(short)*len(pars))
    assert(fc.read!=NULL)

    # fill parameter file
    i = 0
    for kk in sorted(pars):

        dumcp = kk.encode()
        strncpy(fc.name[i], dumcp[:sizeof(FileArg)-1], sizeof(FileArg))
        dumcp = val2str(pars[kk]).encode()
        strncpy(fc.value[i], dumcp[:sizeof(FileArg)-1], sizeof(FileArg))
        fc.read[i] = _FALSE_

        i += 1

    return 0

cdef np.dtype _titles_to_dtype(char * titles, int remove_units=False):
    tmp = (<bytes>titles).decode()
    names = tmp.split("\t")[:-1]
    number_of_titles = len(names)
    if remove_units:
        dtype = np.dtype([(str(name.split()[0]), 'f8') for name in names])
    else:
        dtype = np.dtype([(str(name), 'f8') for name in names])
    return dtype



def _build_task_dependency(tasks):
    """
    Fill the tasks list with all the needed modules

    .. warning::

        the ordering of modules is obviously dependent on CLASS module order
        in the main.c file. This has to be updated in case of a change to
        this file.

    Parameters
    ----------

    tasks : list
        list of strings, containing initially only the last module required.
        For instance, to recover all the modules, the input should be
        ['lensing']

    """

    if "lensing" in tasks:
        tasks.append("spectra")
    if "spectra" in tasks:
        tasks.append("transfer")
    if "transfer" in tasks:
        tasks.append("nonlinear")
    if "nonlinear" in tasks:
        tasks.append("primordial")
    if "primordial" in tasks:
        tasks.append("perturb")
    if "perturb" in tasks:
        tasks.append("thermodynamics")
    if "thermodynamics" in tasks:
        tasks.append("background")
    if len(tasks)!=0 :
        tasks.append("input")
    return tasks


ctypedef struct ready_flags:
    int fc
    int ba
    int th
    int pt
    int pm
    int nl
    int tr
    int sp
    int op
    int le
    int input


cdef class ClassEngine:

    cdef precision pr
    cdef background ba
    cdef thermo th
    cdef perturbs pt
    cdef primordial pm
    cdef nonlinear nl
    cdef transfers tr
    cdef spectra sp
    cdef output op
    cdef lensing le
    cdef ready_flags ready
    cdef file_content fc

    property parameter_file:
        def __get__(self):
            if not self.ready.fc: return ""

            lines = ["%s : %s" %(self.fc.name[i].decode(), self.fc.value[i].decode())
               for i in range(self.fc.size)]
            return "\n".join(lines)

    def __cinit__(self, *args, **kwargs):
        memset(&self.ready, 0, sizeof(self.ready))

    def __init__(self, object pars={}):
        _build_file_content(pars, &self.fc)
        self.ready.fc = True
        self.compute('input')

    def __dealloc__(self):
        if self.ready.fc: parser_free(&self.fc)
        if self.ready.le: lensing_free(&self.le)
        if self.ready.sp: spectra_free(&self.sp)
        if self.ready.tr: transfer_free(&self.tr)
        if self.ready.nl: nonlinear_free(&self.nl)
        if self.ready.pm: primordial_free(&self.pm)
        if self.ready.pt: perturb_free(&self.pt)
        if self.ready.th: thermodynamics_free(&self.th)
        if self.ready.ba: background_free(&self.ba)

    cdef compute(self, level):
        """
        Main function, execute all the _init methods for all desired modules.
        This is called in MontePython, and this ensures that the Class instance
        of this class contains all the relevant quantities. Then, one can deduce
        Pk, Cl, etc...

        Parameters
        ----------
        level : level of modules to arrive.


        """
        cdef file_content * fc = &self.fc
        cdef ErrorMsg errmsg

        tasks = _build_task_dependency([level])

        # --------------------------------------------------------------------
        # Check the presence for all CLASS modules in the list 'tasks'. If a
        # module is found in tasks, executure its "_init" method.
        # --------------------------------------------------------------------
        # The input module should raise a ClassRuntimeError, because
        # non-understood parameters asked to the wrapper is a problematic
        # situation.
        if "input" in tasks and not self.ready.input:
            if input_init(fc, &self.pr, &self.ba, &self.th,
                          &self.pt, &self.tr, &self.pm, &self.sp,
                          &self.nl, &self.le, &self.op, errmsg) == _FAILURE_:
                raise ClassRuntimeError(errmsg.decode())

            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            for i in range(fc.size):
                if fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(fc.name[i].decode())

            if problem_flag:
                import warnings
                warnings.warn("Class did not read input parameter(s): %s" % ', '.join(
                              problematic_parameters))
            self.ready.input = True

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a ClassBadValueError
        # with the error message from the faulty module of CLASS.
        if "background" in tasks and not self.ready.ba:
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                raise ClassBadValueError(self.ba.error_message.decode())
            self.ready.ba = True

        if "thermodynamics" in tasks and not self.ready.th:
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                raise ClassBadValueError(self.th.error_message.decode())
            self.ready.th = True

        if "perturb" in tasks and not self.ready.pt:
            if perturb_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                raise ClassBadValueError(self.pt.error_message.decode())
            self.ready.pt = True

        if "primordial" in tasks and not self.ready.pm:
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                raise ClassBadValueError(self.pm.error_message.decode())
            self.ready.pm = True

        if "nonlinear" in tasks and not self.ready.nl:
            if nonlinear_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.nl) == _FAILURE_:
                raise ClassBadValueError(self.nl.error_message.decode())
            self.ready.nl = True

        if "transfer" in tasks and not self.ready.tr:
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.nl), &(self.tr)) == _FAILURE_:
                raise ClassBadValueError(self.tr.error_message.decode())
            self.ready.tr = True

        if "spectra" in tasks and not self.ready.sp:
            if spectra_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.nl), &(self.tr),
                            &(self.sp)) == _FAILURE_:
                raise ClassBadValueError(self.sp.error_message.decode())
            self.ready.sp = True

        if "lensing" in tasks and not self.ready.le:
            if lensing_init(&(self.pr), &(self.pt), &(self.sp),
                            &(self.nl), &(self.le)) == _FAILURE_:
                raise ClassBadValueError(self.le.error_message.decode())
            self.ready.le = True

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return

    @classmethod
    def from_astropy(cls, cosmo, extra={}):
        """
        Convert an astropy cosmology to a ``ClassEngine`` instance.

        .. warning::

            This method is deprecated. Use classylss.cosmology.Cosmology instead.
        """
        import warnings
        warnings.warn("creating ClassEngine directly from Astropy Cosmology is deprecated; Use classylss.cosmology.Cosmology", DeprecationWarning)

        from astropy import units, cosmology

        pars = {}
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
            if hasattr(m_nu, 'unit') and m_nu.unit != units.eV:
                m_nu = m_nu.to(units.eV)
            else:
                m_nu = units.eV * m_nu
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

        # handle dark energy
        if isinstance(cosmo, cosmology.FlatLambdaCDM):
            #pars['w0_fld'] = -1
            #pars['wa_fld'] = 0.
            #pars['Omega_Lambda'] = 0
            pass
        elif isinstance(cosmo, cosmology.FlatwCDM):
            pars['w0_fld'] = cosmo.w0
            pars['wa_fld'] = 0.
            pars['Omega_Lambda'] = 0. # use Omega_fld
        elif isinstance(cosmo, cosmology.Flatw0waCDM):
            pars['w0_fld'] = cosmo.w0
            pars['wa_fld'] = cosmo.wa
            pars['Omega_Lambda'] = 0. # use Omega_fld
        else:
            cls = cosmo.__class__.__name__
            valid = ["FlatLambdaCDM", "FlatwCDM", "Flatw0waCDM"]
            msg = "dark energy equation of state not recognized for class '%s'; " %cls
            msg += "valid classes: %s" %str(valid)
            raise ValueError(msg)

        # add any extra arguments
        if len(extra):
            pars.update(extra)

        _pars = {
            "output": "vTk dTk mPk",
            "P_k_max_h/Mpc":  100.,
            "z_max_pk": 100.0,
            "extra metric transfer functions": "y",
            }

        for k in _pars:
            pars.setdefault(k, _pars[k])

        self = cls(pars)
        return self

cdef class Background:
    cdef ClassEngine engine
    cdef background * ba
    cdef readonly dict data

    cdef readonly double Omega0_pncdm
    cdef readonly double H0
    cdef readonly double C
    cdef readonly double G
    cdef readonly double _RHO_

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("background")
        self.ba = &self.engine.ba

        self.H0 = 100.  # in Mpc/h unit
        self.G = 43007.1 * 1e-3 # in 1e10 Msun/h, Mpc/h, and km/s Unit
        self.C = 2.99792458e5          #  /**< c in km/s */
        # convert RHO to  1e10 Msun/h
        self._RHO_ = 3.0 * (self.H0 / self.ba.H0) ** 2 / (8 * 3.1415927 * self.G)

        self.Omega0_pncdm = self.Omega_pncdm(0.0) # watchout, the convention is 0.0

    property Omega0_b:
        """
        Density parameter for photons, :math:`\Omega_{b,0}`
        """
        def __get__(self):
            return self.ba.Omega0_b

    property Omega0_g:
        """
        Density parameter for photons, :math:`\Omega_{g,0}`
        """
        def __get__(self):
            return self.ba.Omega0_g

    property Omega0_cdm:
        """
        Density parameter for cold dark matter, :math:`\Omega_{cdm,0}`
        """
        def __get__(self):
            return self.ba.Omega0_cdm

    property Omega0_lambda:
        """
        Density parameter for cosmological constant, :math:`\Omega_{\Lambda,0}`
        """
        def __get__(self):
            return self.ba.Omega0_lambda

    property Omega0_fld:
        """
        Density parameter for dark energy (fluid), :math:`\Omega_{de,0}`
        """
        def __get__(self):
            return self.ba.Omega0_fld

    property Omega0_k:
        """
        Density parameter for curvaturve, :math:`\Omega_{k,0}`
        """
        def __get__(self):
            return self.ba.Omega0_k

    property w0_fld:
        """
        Current fluid equation of state parameter, :math:`w_{0,fld}`
        """
        def __get__(self):
            return self.ba.w0_fld

    property wa_fld:
        """
        Fluid equation of state derivative, :math:`w_{a,fld}`
        """
        def __get__(self):
            return self.ba.wa_fld

    property Omega0_dcdm:
        """
        Density parammeter for decaying dark matter, :math:`\Oemga_{dcdm,0}`
        """
        def __get__(self):
            return self.ba.Omega0_dcdm

    property Omega0_ncdm:
        """
        Density of distinguishable (massive) neutrinos for each species as an array
        """
        def __get__(self):
            return np.array([self.ba.Omega0_ncdm[i] for i in range(self.N_ncdm)], np.float64)

    property Omega0_ncdm_tot:
        """
        Total density of distinguishable (massive) neutrinos ;
        """
        def __get__(self):
            return self.ba.Omega0_ncdm_tot

    property Omega0_ur:
        """
        Density parameter of ultra-relativistic (massless) neutrinos,
        :math:`\Omega_{0,\nu_r}`
        """
        def __get__(self):
            return self.ba.Omega0_ur

    property Omega0_r:
        def __get__(self):
            return self.ba.Omega0_g + self.ba.Omega0_ur + self.Omega0_pncdm

    property a_today:
        """ this is an arbitrary number that sets the refernce scaling factor. It shall be 1 usually."""
        def __get__(self):
            return self.ba.a_today

    property a_max:
        def __get__(self):
            return self.ba.a_max

    property Omega0_m:
        """
        Return the sum of Omega0 for all non-relativistic components. The value differ from Astropy's; the semantics is identical.
        """
        def __get__(self):
            return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm - self.Omega0_pncdm

    property Neff:
        def __get__(self):
            return self.ba.Neff

    property N_ur:
        def __get__(self):
            return self.Omega0_ur / (7./8.*(4./11)**(4./3.)*self.Omega0_g)

    property N_ncdm:
        def __get__(self):
            return self.ba.N_ncdm

    property m_ncdm:
        def __get__(self):
            return np.array([self.ba.m_ncdm_in_eV[i] for i in range(self.N_ncdm)], dtype=np.float64)

    property age0:
        def __get__(self):
            return self.ba.age

    property h:
        def __get__(self):
            return self.ba.h

    property T0_cmb:
        """
        Return the CMB temperature
        """
        def __get__(self):
            return self.ba.T_cmb

    def compute_for_z(self, z, int column):
        cdef double tau
        cdef int last_index #junk

        z = np.array(z, dtype=np.float64)

        #generate a new output array of the correct shape by broadcasting input arrays together
        out = np.empty(np.broadcast(z).shape, np.float64)

        #generate the iterator over the input and output arrays, does the same thing as
        # PyArray_MultiIterNew

        cdef np.broadcast it = np.broadcast(z, out)

        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')

        while np.PyArray_MultiIter_NOTDONE(it):

            #PyArray_MultiIter_DATA is used to access the pointers the iterator points to
            aval = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]

            if background_tau_of_z(self.ba,aval, &tau)==_FAILURE_:
                raise ClassRuntimeError(self.ba.error_message.decode())

            if background_at_tau(self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index, &pvecback[0])==_FAILURE_:
                raise ClassRuntimeError(self.ba.error_message.decode())

            (<double*>(np.PyArray_MultiIter_DATA(it, 1)))[0] = pvecback[column]

            #PyArray_MultiIter_NEXT is used to advance the iterator
            np.PyArray_MultiIter_NEXT(it)

        return out

    def Omega_pncdm(self, z, species=None):
        return 3 * self.p_ncdm(z, species) / self.rho_tot(z)

    def rho_g(self, z):
        """ density of radiation (gamma). """
        return self.compute_for_z(z, self.ba.index_bg_rho_g) * self._RHO_

    def rho_b(self, z):
        """ density of baryon. """
        return self.compute_for_z(z, self.ba.index_bg_rho_b) * self._RHO_

    def rho_m(self, z):
        """ density of matter like components."""
        return self.Omega_m(z) * self.rho_tot(z)

    def rho_r(self, z):
        """ density of radiation like components. """
        return self.Omega_r(z) * self.rho_tot(z)

    def rho_cdm(self, z):
        """ density of cdm. """
        return self.compute_for_z(z, self.ba.index_bg_rho_cdm) * self._RHO_

    def rho_ur(self, z):
        """ density of ultra relative (massless) neutrino. """
        return self.compute_for_z(z, self.ba.index_bg_rho_ur) * self._RHO_

    def rho_ncdm(self, z, species=None):
        """ density of non-relative part of massive neutrino. """
        if species is None:
            return sum(self.rho_ncdm(z, species=i) for i in range(self.N_ncdm))
        assert species < self.N_ncdm and species >= 0
        return self.compute_for_z(z, self.ba.index_bg_rho_ncdm1 + species) * self._RHO_

    def rho_crit(self, z):
        """ critical density, total excluding curvature """
        return self.compute_for_z(z, self.ba.index_bg_rho_crit) * self._RHO_

    def rho_k(self, z):
        """ curvature density. """
        z = np.array(z, dtype=np.float64)
        return self.ba.K * ( z+1.) ** 2 * self._RHO_

    def rho_tot(self, z):
        """ total density in 1e10 Msun/h (Mpc/h) ** -3; it is usually close to 27.76."""
        return self.rho_crit(z) + self.rho_k(z)

    def rho_fld(self, z):
        """ density of fluid (dark energy fluid). """
        if self.ba.has_fld:
            return self.compute_for_z(z, self.ba.index_bg_rho_fld) * self._RHO_
        else:
            # return zeros of the right shape
            return self.compute_for_z(z, self.ba.index_bg_a) * 0.0

    def rho_lambda(self, z):
        """ density of cosmology constant (dark energy non fluid) """
        if self.ba.has_lambda:
            return self.compute_for_z(z, self.ba.index_bg_rho_lambda) * self._RHO_
        else:
            # return zeros of the right shape
            return self.compute_for_z(z, self.ba.index_bg_a) * 0.0

    def p_ncdm(self, z, species=None):
        """ pressure of non-relative part of massive neutrino. """
        if species is None:
            return sum(self.p_ncdm(z, i) for i in range(self.N_ncdm))

        assert species < self.N_ncdm and species >= 0
        return self.compute_for_z(z, self.ba.index_bg_p_ncdm1 + species) * self._RHO_

    def Omega_r(self, z):
        """ density of relative (radiation like) component, including relative part of massive neutrino and massless neutrino. """
        return self.compute_for_z(z, self.ba.index_bg_Omega_r)

    def Omega_m(self, z):
        """ density of non-relative (matter like) component, including non-relative part of massive neutrino. """
        return self.compute_for_z(z, self.ba.index_bg_Omega_m)

    def Omega_g(self, z):
        """ density of radiation """
        return self.rho_g(z) / self.rho_tot(z)

    def Omega_b(self, z):
        """ density of baryon """
        return self.rho_b(z) / self.rho_tot(z)

    def Omega_cdm(self, z):
        """ density of baryon """
        return self.rho_cdm(z) / self.rho_tot(z)

    def Omega_k(self, z):
        """ density of curvature"""
        return self.rho_k(z) / self.rho_tot(z)

    def Omega_ur(self, z):
        """ density of ultra relativistic neutrino """
        return self.rho_ur(z) / self.rho_tot(z)

    def Omega_ncdm(self, z, species=None):
        """ density of massive neutrino """
        return self.rho_ncdm(z, species) / self.rho_tot(z)

    def Omega_fld(self, z):
        """ density of dark energy (fluid) """
        return self.rho_fld(z) / self.rho_tot(z)

    def Omega_lambda(self, z):
        """ density of dark energy (cosmological constant) """
        return self.rho_lambda(z) / self.rho_tot(z)

    def time(self, z):
        """ proper time (age of universe) """
        return self.compute_for_z(z, self.ba.index_bg_time)

    def conformal_distance(self, z):
        """ conformal distance, comoving distance when K = 0.0 (flat universe). In units of Mpc/h"""
        return self.compute_for_z(z, self.ba.index_bg_conf_distance) * self.ba.h

    def tau(self, z):
        """ conformal time , comoving distance when K = 0.0 (flat universe). In units of Mpc as CLASS"""
        return self.compute_for_z(z, self.ba.index_bg_conf_distance)

    def hubble_function(self, z):
        """ hubble function in class units; use efunc instead """
        return self.compute_for_z(z, self.ba.index_bg_H)

    def hubble_function_prime(self, z):
        """ d H / d tau ; d tau / da = 1 / (a ** 2 H) in class units; use efunc_prime instead """
        return self.compute_for_z(z, self.ba.index_bg_H_prime)

    def efunc(self, z):
        """
        Function giving :func:`efunc with respect
        to the scale factor ``a``
        Parameters
        ----------
        z : array-like
            Input redshifts.
        Returns
        -------
        efunc : ndarray, or float if input scalar
            The hubble factor redshift-scaling with respect
            to scale factor
        """
        return self.hubble_function(z) / self.ba.H0

    def efunc_prime(self, z):
        """
        Function giving dE / da.
        """
        dtau_da = (1 + z)**2 / self.hubble_function(z)
        return self.hubble_function_prime(z) / self.ba.H0 * dtau_da

    def luminosity_distance(self, z):
        """ in Mpc/h"""
        return self.compute_for_z(z, self.ba.index_bg_lum_distance) * self.ba.h

    def angular_distance(self, z):
        """
        in Mpc/h

        Parameters
        ----------
        z : float
                Desired redshift
        """
        return self.compute_for_z(z, self.ba.index_bg_ang_distance) * self.ba.h

    def scale_independent_growth_factor(self, z):
        """
        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        return self.compute_for_z(z, self.ba.index_bg_D)

    def scale_independent_growth_rate(self, z):
        """

        Return the scale invariant growth rate d ln D(a) / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        return self.compute_for_z(z, self.ba.index_bg_f)

cdef class Perturbs:
    cdef ClassEngine engine
    cdef perturbs * pt
    cdef background * ba

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("perturbs")
        self.pt = &self.engine.pt
        self.ba = &self.engine.ba

    property P_k_max:
        def __get__(self):
            return self.pt.k_max_for_pk/self.ba.h

    property P_z_max:
        def __get__(self):
            return self.pt.z_max_pk

    property gauge:
        def __get__(self):
            return self.pt.gauge

cdef class Primordial:
    cdef ClassEngine engine
    cdef perturbs * pt
    cdef primordial * pm
    cdef background * ba

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("primordial")
        self.pt = &self.engine.pt
        self.ba = &self.engine.ba
        self.pm = &self.engine.pm

    def get_pk(self, k, mode='linear'):
        """ get primoridal spectrum at k.

            Parameters
            ----------
            k : array_like,  h/Mpc units.

            Results
            -------
            primordial_k : array_like

                Unit unclear. Tf_k ** 2 * primordial_k is dimensionless.

        """

        #generate a new output array of the correct shape by broadcasting input arrays together
        k = np.float64(k) * self.ba.h
        out = np.empty(np.broadcast(k).shape, np.float64)

        #generate the iterator over the input and output arrays, does the same thing as
        cdef np.broadcast it = np.broadcast(k,  out)
        cdef int index_md = 0
        cdef linear_or_logarithmic modeval

        if mode.startswith('lin'):
            modeval = linear
        elif mode.startswith('log'):
            modval = logarithmic
        else:
            raise ValueError("mode can only be log or lin")

        while np.PyArray_MultiIter_NOTDONE(it):

                #PyArray_MultiIter_DATA is used to access the pointers the iterator points to
                aval = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]

                if _FAILURE_ == primordial_spectrum_at_k(self.pm, index_md, modeval, aval,
                    <double*>(np.PyArray_MultiIter_DATA(it, 1))):
                    raise ClassRuntimeError(self.pm.error_message.decode())

                #PyArray_MultiIter_NEXT is used to advance the iterator
                np.PyArray_MultiIter_NEXT(it)

        # Watch out: no transformation here
        return out

    def get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """

        primordial = {}
        cdef char titles[_MAXTITLESTRINGLENGTH_]
        memset(titles, 0, _MAXTITLESTRINGLENGTH_)

        if primordial_output_titles(self.pt, self.pm, titles)==_FAILURE_:
            raise ClassRuntimeError(self.pm.error_message.decode())

        dtype = _titles_to_dtype(titles)

        cdef np.ndarray data = np.zeros(self.pm.lnk_size, dtype=dtype)

        if primordial_output_data(self.pt, self.pm, len(dtype.fields), <double*>data.data)==_FAILURE_:
            raise ClassRuntimeError(self.pm.error_message.decode())

        return data

cdef class Spectra:
    cdef ClassEngine engine
    cdef spectra * sp
    cdef background * ba
    cdef perturbs * pt
    cdef primordial * pm
    cdef nonlinear * nl
    cdef readonly dict data

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("spectra")

        self.ba = &self.engine.ba
        self.nl = &self.engine.nl
        self.sp = &self.engine.sp
        self.pt = &self.engine.pt
        self.pm = &self.engine.pm

    property nonlinear:
      def __get__(self):
        return self.nl.method > 0

    property has_pk_matter:
        def __get__(self):
          return self.pt.has_pk_matter

    property P_k_min:
        def __get__(self):
            return np.exp(self.sp.ln_k[0])/self.ba.h;

    property sigma8:
        def __get__(self):
            return self.sp.sigma8

    property A_s:
        def __get__(self):
            return self.pm.A_s

    property ln_1e10_A_s:
        def __get__(self):
            return np.log(1e10*self.A_s)

    property n_s:
        def __get__(self):
            return self.pm.n_s

    def sigma8_z(self, z):

      #generate a new output array of the correct shape by broadcasting input arrays together
      z = np.float64(z)
      out = np.empty(np.broadcast(z).shape, np.float64)

      #generate the iterator over the input and output arrays, does the same thing as
      cdef np.broadcast it = np.broadcast(z,  out)

      while np.PyArray_MultiIter_NOTDONE(it):

              #PyArray_MultiIter_DATA is used to access the pointers the iterator points to
              aval = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]

              if _FAILURE_ == spectra_sigma(self.ba, self.pm, self.sp, 8./self.ba.h, aval,
                  <double*>(np.PyArray_MultiIter_DATA(it, 1))):
                  raise ClassRuntimeError(self.sp.error_message.decode())

              #PyArray_MultiIter_NEXT is used to advance the iterator
              np.PyArray_MultiIter_NEXT(it)

      return out

    def get_transfer(self, z, output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'dCl' and 'vCl' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that z is inside the region spanned by 'z_pk'.

        This function is not vectorized; because it returns a vector when ic_size is
        greater than 1, and I don't understand ic_size.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : array_like, containing transfer functions. Unlike CLASS, k here is in Mpc/h Units.

        .. note::

            With different cosmology the values of 'k' may be different.
        """

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            raise RuntimeError("Perturbation is not computed")

        cdef FileName ic_suffix
        cdef file_format_outf
        cdef char ic_info[1024]

        cdef char titles[_MAXTITLESTRINGLENGTH_]
        memset(titles, 0, _MAXTITLESTRINGLENGTH_)

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        if spectra_output_tk_titles(self.ba, self.pt, outf, titles)==_FAILURE_:
            raise ClassRuntimeError(self.op.error_message.decode())

        # k is in h/Mpc. Other functions unit is unclear.
        dtype = _titles_to_dtype(titles, remove_units=True)

        index_md = 0
        ic_num = self.sp.ic_size[index_md]

        cdef np.ndarray data = np.zeros((ic_num, self.sp.ln_k_size), dtype=dtype)

        if spectra_output_tk_data(self.ba, self.pt, self.sp, outf, <double> z, len(dtype.fields), <double*> data.data)==_FAILURE_:
            raise ClassRuntimeError(self.sp.error_message.decode())

        ic_keys = []
        if ic_num > 1:
            for index_ic in range(ic_num):
                if spectra_firstline_and_ic_suffix(self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                    raise ClassRuntimeError(self.op.error_message.decode())

                ic_key = <bytes> ic_suffix
                ic_keys.append(ic_key.decode())

            spectra = {}
            for ic_key, row in zip(ic_keys, data):
                spectra[ic_key] = row
        else:
            spectra = data[0]

        return spectra


    # Gives the pk for a given (k,z)
    cdef int pk(self, double k, double z, double * pk_ic, int lin, double * pk) except -1:
        """
        Gives the pk for a given k and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        if lin or self.nl.method == 0:
            if spectra_pk_at_k_and_z(self.ba,self.pm,self.sp,k,z,pk,pk_ic) == _FAILURE_:
                 raise ClassRuntimeError(self.sp.error_message.decode())
        else:
            if spectra_pk_nl_at_k_and_z(self.ba,self.pm,self.sp,k,z,pk) ==_FAILURE_:
                 raise ClassRuntimeError(self.sp.error_message.decode())
        return 0

    def get_pk(self, k, z):
        """ Primary Power spectrum result (non-linear if enabled) on k and z array. K in h/Mpc units.

            Results
            -------
            Pk : array like
                Power Spectrum in (Mpc/h)**3
        """
        return self._get_pk(k, z, 0)

    def get_pklin(self, k, z):
        """ Linear Power spectrum result (linear even if non-linear is enabled) on k and z array. K in h/Mpc units.

            Results
            -------
            Pk : array like
                Power Spectrum in (Mpc/h)**3
        """
        return self._get_pk(k, z, 1)

    def _get_pk(self, k, z, int linear):

        if (self.pt.has_pk_matter == _FALSE_):
            raise ClassRuntimeError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )

        # internally class uses 1 / Mpc
        k = np.float64(k) * self.ba.h
        z = np.float64(z)

        # Quantities for the isocurvature modes
        cdef np.ndarray pk_ic = np.zeros(self.sp.ic_ic_size[self.sp.index_md_scalars], dtype='f8')

        #generate a new output array of the correct shape by broadcasting input arrays together
        out = np.empty(np.broadcast(k, z).shape, np.float64)

        #generate the iterator over the input and output arrays, does the same thing as
        # PyArray_MultiIterNew

        cdef np.broadcast it = np.broadcast(k, z, out)

        while np.PyArray_MultiIter_NOTDONE(it):

                #PyArray_MultiIter_DATA is used to access the pointers the iterator points to
                aval = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]
                bval = (<double*>np.PyArray_MultiIter_DATA(it, 1))[0]

                self.pk(aval, bval, <double*> pk_ic.data, linear, <double*>(np.PyArray_MultiIter_DATA(it, 2)))

                #PyArray_MultiIter_NEXT is used to advance the iterator
                np.PyArray_MultiIter_NEXT(it)

        # internally class uses Mpc ** 3
        out[...] *= self.ba.h**3

        return out
