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

class ClassParserError(ValueError):
    def __init__(self, message, file_content):
        self.message = message
        self.file_content = file_content

    def __str__(self):
        return 'Class Parser Error : ' + self.message + '\n' + self.file_content

    pass

class ClassBadValueError(ValueError):
    r"""
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass


def val2str(val):
    if isinstance(val, (list, tuple, np.ndarray)):
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
    r"""
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
    """
    The default CLASS engine class, which initializes CLASS from an input
    set of parameters.

    Parameters
    ----------
    pars : dict, optional
      a dictionary of parameters to initialize CLASS with
    """
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
        """
        A string holding the parameter names and values as loaded by CLASS.
        """
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
        r"""
        The main function, which executes all the 'init' methods for all
        the desired modules.

        Parameters
        ----------
        level : str
          level of modules to arrive.
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
                raise ClassParserError(errmsg.decode(), self.parameter_file)

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

cdef class Background:
    """
    A wrapper of the `background module <https://goo.gl/SU71dn>`_ in CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
    cdef ClassEngine engine
    cdef background * ba
    cdef readonly dict data

    cdef readonly np.ndarray Omega0_pncdm
    """
    The pressure contribution to the current density parameter for the
    non-relativatistic part of massive neutrinos (an array holding all species).
    """
    cdef readonly double Omega0_pncdm_tot
    """
    The sum of :math:`\Omega_{0,pncdm}` for all species.
    """
    cdef readonly double H0
    r"""
    Current Hubble parameter in units of :math:`\mathrm{km/s} (\mathrm{Mpc}/h)^{-1}.`
    """
    cdef readonly double C
    """
    The speed of light in units of km/s.
    """
    cdef readonly double G
    r"""
    The gravitational constant in units of :math:`10^{-10} \ (M_\odot/h)^{-1} (\mathrm{Mpc}/h) \mathrm{km}^2 \mathrm{s}^{-2}`.
    """

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

        self.Omega0_pncdm_tot = self.Omega_pncdm(0.0) # watchout, the convention is 0.0
        self.Omega0_pncdm = np.array([self.Omega_pncdm(0.0, i) for i in range(self.N_ncdm)], np.float64)

    property Omega0_b:
        r"""
        Current density parameter for photons, :math:`\Omega_{b,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_b

    property Omega0_g:
        r"""
        Current density parameter for photons, :math:`\Omega_{g,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_g

    property Omega0_cdm:
        r"""
        Current density parameter for cold dark matter, :math:`\Omega_{cdm,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_cdm

    property Omega0_lambda:
        r"""
        Current density parameter for cosmological constant, :math:`\Omega_{\Lambda,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_lambda

    property Omega0_fld:
        r"""
        Current density parameter for dark energy (fluid) :math:`\Omega_{fld, 0}`.
        """
        def __get__(self):
            return self.ba.Omega0_fld

    property Omega0_k:
        r"""
        Current density parameter for curvaturve, :math:`\Omega_{k,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_k

    property w0_fld:
        r"""
        Current fluid equation of state parameter, :math:`w_{0,fld}`.
        """
        def __get__(self):
            return self.ba.w0_fld

    property wa_fld:
        r"""
        Fluid equation of state derivative, :math:`w_{a,fld}`.
        """
        def __get__(self):
            return self.ba.wa_fld

    property Omega0_dcdm:
        r"""
        Current density parammeter for decaying cold dark matter,
        :math:`\Omega_{dcdm,0}`.
        """
        def __get__(self):
            return self.ba.Omega0_dcdm

    property Omega0_ncdm:
        r"""
        Current density parameter for distinguishable (massive) neutrinos for
        each species as an array, :math:`\Omega_{0, ncdm}`.
        """
        def __get__(self):
            return np.array([self.ba.Omega0_ncdm[i] for i in range(self.N_ncdm)], np.float64)

    property Omega0_ncdm_tot:
        r"""
        Current total density parameter of all distinguishable (massive)
        neutrinos.
        """
        def __get__(self):
            return self.ba.Omega0_ncdm_tot

    property Omega0_ur:
        r"""
        Current density parameter of ultra-relativistic (massless) neutrinos,
        :math:`\Omega_{0,\nu_r}`.
        """
        def __get__(self):
            return self.ba.Omega0_ur

    property Omega0_r:
        r"""
        Current density parameter of radiation, :math:`\Omega_{0,r}`.
        This is equal to:

        .. math::

            \Omega_{0,r} = \Omega_{0,g} + \Omega_{0,\nu_r} + \Omega_{0,pncdm}.
        """
        def __get__(self):
            return self.ba.Omega0_g + self.ba.Omega0_ur + self.Omega0_pncdm_tot

    property a_today:
        r"""
        An arbitrary number that sets the reference scaling factor.
        It shall be 1 usually.
        """
        def __get__(self):
            return self.ba.a_today

    property a_max:
        r"""
        The maximum scale factor for which results can be computed; it can be
        greater than 1.0.
        """
        def __get__(self):
            return self.ba.a_max

    property Omega0_m:
        r"""
        The sum of density parameters for all non-relativistic components,
        :math:`\Omega_{0,m}`. The value differs from ``Om0`` in :mod:`astropy`.

        This is equal to:

        .. math::
            \Omega_{0,m} = \Omega_{0,b} + \Omega_{0,cdm} + \Omega_{0,ncdm} + \Omega_{0,dcdm} - \Omega_{0,pncdm}.
        """
        def __get__(self):
            return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + \
                  self.ba.Omega0_dcdm - self.Omega0_pncdm_tot

    property Neff:
        r"""
        Effective number of relativistic species, summed over ultra-relativistic
        and ncdm species.
        """
        def __get__(self):
            return self.ba.Neff

    property N_ur:
        r"""
        The number of ultra-relativistic species.

        This is equal to:

        .. math::

            N_{ur} = \Omega_{0,ur} / (7/8 (4/11)^{4/3} \Omega_{0,g}).
        """
        def __get__(self):
            return self.Omega0_ur / (7./8.*(4./11)**(4./3.)*self.Omega0_g)

    property N_ncdm:
        r"""
        The number of distinguishable ncdm (massive neutrino) species.
        """
        def __get__(self):
            return self.ba.N_ncdm

    property m_ncdm:
        r"""
        The masses of the distinguishable ncdm (massive neutrino) species,
        in units of eV.
        """
        def __get__(self):
            return np.array([self.ba.m_ncdm_in_eV[i] for i in range(self.N_ncdm)], dtype=np.float64)

    property age0:
        r"""
        The current age of the universe in gigayears.
        """
        def __get__(self):
            return self.ba.age

    property h:
        r"""
        The dimensionless Hubble parameter.
        """
        def __get__(self):
            return self.ba.h

    property T0_cmb:
        r"""
        The current CMB temperature in Kelvins.
        """
        def __get__(self):
            return self.ba.T_cmb

    property T0_ncdm:
        r"""
        An array holding the current ncdm temperature in Kelvins for each species.
        """
        def __get__(self):
            T = np.array([self.ba.T_ncdm[i] for i in range(self.N_ncdm)], dtype=np.float64)
            return T*self.ba.T_cmb # from units of photon temp to K

    def T_cmb(self, z):
        r"""
        The CMB temperature as a function of redshift.
        """
        return self.T0_cmb * (1 + z)

    def T_ncdm(self, z):
        r"""
        The ncdm temperature (massive neutrinos) as a function of redshift.

        Return shape is (N_ncdm,) if N_ncdm == 1 else (len(z), N_ncdm)
        """
        if np.isscalar(z):
            return self.T0_ncdm * (1+z)
        else:
            z = np.array(z, ndmin=1, dtype=np.float64)
        return self.T0_ncdm * (1 + z)[:,None]

    def compute_for_z(self, z, int column):
        """
        Internal function to compute the background module at a specific redshift.
        """
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
        r"""
        Return :math:`\Omega_{pncdm}` as a function redshift.
        """
        return 3 * self.p_ncdm(z, species) / self.rho_tot(z)

    def rho_g(self, z):
        r"""
        Density of photons :math:`\rho_g` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_rho_g) * self._RHO_

    def rho_b(self, z):
        r"""
        Density of baryons :math:`\rho_b` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_rho_b) * self._RHO_

    def rho_m(self, z):
        r"""
        Density of matter :math:`\rho_b` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_Omega_m) * self.rho_tot(z)

    def rho_r(self, z):
        r"""
        Density of radiation :math:`\rho_r` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_Omega_r) * self.rho_tot(z)

    def rho_cdm(self, z):
        r"""
        Density of cold dark matter :math:`\rho_{cdm}` as a function of redshift,
        in units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_rho_cdm) * self._RHO_

    def rho_ur(self, z):
        r"""
        Density of ultra-relativistic radiation (massless neutrinos)
        :math:`\rho_{ur}` as a function of redshift, in units of
        :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        return self.compute_for_z(z, self.ba.index_bg_rho_ur) * self._RHO_

    def rho_ncdm(self, z, species=None):
        r"""
        Density of non-relativistic part of massive neutrinos :math:`\rho_{ncdm}`
        as a function of redshift, in units of
        :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        if species is None:
            return sum(self.rho_ncdm(z, species=i) for i in range(self.N_ncdm))
        assert species < self.N_ncdm and species >= 0
        return self.compute_for_z(z, self.ba.index_bg_rho_ncdm1 + species) * self._RHO_

    def rho_crit(self, z):
        r"""
        Critical density excluding curvature :math:`\rho_c` as a function of
        redshift, in units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.

        This is defined as:

        .. math::

              \rho_c(z) = \frac{3 H(z)^2}{8 \pi G}.
        """
        return self.compute_for_z(z, self.ba.index_bg_rho_crit) * self._RHO_

    def rho_k(self, z):
        r"""
        Density of curvature :math:`\rho_k` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.

        Note: this is defined such that

        .. math::

            \rho_\mathrm{crit} = \rho_\mathrm{tot} + \rho_k
        """
        z = np.array(z, dtype=np.float64)
        return -self.ba.K * ( z+1.) ** 2 * self._RHO_

    def rho_tot(self, z):
        r"""
        Total density :math:`\rho_\mathrm{tot}` as a function of redshift, in
        units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`. It is usually
        close to 27.76.
        """
        return self.rho_crit(z) - self.rho_k(z)

    def rho_fld(self, z):
        r"""
        Density of dark energy fluid :math:`\rho_{fld}` as a function of
        redshift, in units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        if self.ba.has_fld:
            return self.compute_for_z(z, self.ba.index_bg_rho_fld) * self._RHO_
        else:
            # return zeros of the right shape
            return self.compute_for_z(z, self.ba.index_bg_a) * 0.0

    def rho_lambda(self, z):
        r"""
        Density of cosmological constant :math:`\rho_\Lambda` as a function of
        redshift, in units of :math:`10^{10} (M_\odot/h) (\mathrm{Mpc}/h)^{-3}`.
        """
        if self.ba.has_lambda:
            return self.compute_for_z(z, self.ba.index_bg_rho_lambda) * self._RHO_
        else:
            # return zeros of the right shape
            return self.compute_for_z(z, self.ba.index_bg_a) * 0.0

    def p_ncdm(self, z, species=None):
        r"""
        Pressure of non-relative part of massive neutrino.
        """
        if species is None:
            return sum(self.p_ncdm(z, i) for i in range(self.N_ncdm))

        assert species < self.N_ncdm and species >= 0
        return self.compute_for_z(z, self.ba.index_bg_p_ncdm1 + species) * self._RHO_

    def Omega_r(self, z):
        r"""
        Density parameter of relativistic (radiation-like) component, including
        relativistic part of massive neutrino and massless neutrino.
        """
        return self.rho_r(z) / self.rho_crit(z)

    def Omega_m(self, z):
        r"""
        Density parameter of non-relativistic (matter-like) component, including
        non-relativistic part of massive neutrino. Unit
        """
        return self.rho_m(z) / self.rho_crit(z)

    def Omega_g(self, z):
        r"""
        Density parameter of photons.
        """
        return self.rho_g(z) / self.rho_crit(z)

    def Omega_b(self, z):
        r"""
        Density parameter of baryons.
        """
        return self.rho_b(z) / self.rho_crit(z)

    def Omega_cdm(self, z):
        r"""
        Density parameter of cold dark matter.
        """
        return self.rho_cdm(z) / self.rho_crit(z)

    def Omega_k(self, z):
        r"""
        Density parameter of curvature.
        """
        return 1 - self.rho_tot(z)/self.rho_crit(z)

    def Omega_ur(self, z):
        r"""
        Density parameter of ultra relativistic neutrinos.
        """
        return self.rho_ur(z) / self.rho_crit(z)

    def Omega_ncdm(self, z, species=None):
        r"""
        Density parameter of massive neutrinos.
        """
        return self.rho_ncdm(z, species) / self.rho_crit(z)

    def Omega_fld(self, z):
        r"""
        Density parameter of dark energy (fluid).
        """
        return self.rho_fld(z) / self.rho_crit(z)

    def Omega_lambda(self, z):
        r"""
        Density of dark energy (cosmological constant).
        """
        return self.rho_lambda(z) / self.rho_crit(z)

    def time(self, z):
        r"""
        Proper time (age of universe) in gigayears.
        """
        return self.compute_for_z(z, self.ba.index_bg_time) / _Gyr_over_Mpc_

    def comoving_distance(self, z):
        r"""
        Comoving line-of-sight distance in :math:`\mathrm{Mpc}/h` at a given
        redshift.

        The comoving distance along the line-of-sight between two
        objects remains constant with time for objects in the Hubble
        flow.

        See eq. 15 of `astro-ph/9905116 <https://arxiv.org/abs/astro-ph/9905116>`_
        for :math:`D_C(z)`.
        """
        return self.compute_for_z(z, self.ba.index_bg_conf_distance) * self.ba.h

    def tau(self, z):
        r"""
        Conformal time, equal to comoving distance when K = 0.0
        (flat universe). In units of :math:`\mathrm{Mpc}` as in CLASS.
        """
        return self.compute_for_z(z, self.ba.index_bg_conf_distance)

    def hubble_function(self, z):
        r"""
        The Hubble function in CLASS units, returning ``ba.index_bg_H``.

        Users should use :func:`efunc` instead.
        """
        return self.compute_for_z(z, self.ba.index_bg_H)

    def hubble_function_prime(self, z):
        r"""
        Derivative of Hubble function: :math:`dH/d\tau`, where
        :math:`d\tau/da = 1 / (a^2 H)` in CLASS units.

        Users should use :func:`efunc_prime` instead.
        """
        return self.compute_for_z(z, self.ba.index_bg_H_prime)

    def efunc(self, z):
        r"""
        Function giving :math:`E(z)`, where the Hubble parameter is defined as
        :math:`H(z) = H_0 E(z)`.
        """
        return self.hubble_function(z) / self.ba.H0

    def efunc_prime(self, z):
        r"""
        Function giving :math:`dE(z) / da`.
        """
        dtau_da = (1 + z)**2 / self.hubble_function(z)
        return self.hubble_function_prime(z) / self.ba.H0 * dtau_da

    def luminosity_distance(self, z):
        r"""
        Luminosity distance in :math:`\mathrm{Mpc}/h` at redshift ``z``.

        This is the distance to use when converting between the
        bolometric flux from an object at redshift ``z`` and its
        bolometric luminosity.

        It is equal to the comoving transverse distance times :math:`1+z`.

        See eq. 21 of `astro-ph/9905116 <https://arxiv.org/abs/astro-ph/9905116>`_
        for :math:`D_L(z)`.
        """
        return self.compute_for_z(z, self.ba.index_bg_lum_distance) * self.ba.h

    def angular_diameter_distance(self, z):
        r"""
        Angular diameter distance in :math:`\mathrm{Mpc}/h` at a given redshift.

        This gives the proper (sometimes called 'physical') transverse
        distance corresponding to an angle of 1 radian for an object
        at redshift ``z``.

        It is equal to the comoving transverse distance divided by :math:`1+z`.

        See eq. 18 of `astro-ph/9905116 <https://arxiv.org/abs/astro-ph/9905116>`_
        for :math:`D_A(z)`.
        """
        return self.compute_for_z(z, self.ba.index_bg_ang_distance) * self.ba.h

    def comoving_transverse_distance(self, z):
        r"""
        Comoving transverse distance in :math:`\mathrm{Mpc}/h` at a given
        redshift.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is
        the same as the comoving distance in a flat universe.

        See eq. 16 of `astro-ph/9905116 <https://arxiv.org/abs/astro-ph/9905116>`_
        for :math:`D_M(z)`.
        """
        # comoving distance if flat (in Mpc)
        toret = self.comoving_distance(z) / self.ba.h

        # positive curvature
        if (self.ba.sgnK == 1):
            toret = np.sin(np.sqrt(self.ba.K)*toret)/np.sqrt(self.ba.K)

        # negative curvature
        if (self.ba.sgnK == -1):
            toret = np.sinh(np.sqrt(-self.ba.K)*toret)/np.sqrt(-self.ba.K)

        # in Mpc/h
        return toret * self.ba.h

    def scale_independent_growth_factor(self, z):
        r"""
        Return the scale invariant growth factor :math:`D(a)` for CDM
        perturbations.

        This is the quantity defined by CLASS as ``index_bg_D`` in the
        background module.
        """
        return self.compute_for_z(z, self.ba.index_bg_D)

    def scale_independent_growth_rate(self, z):
        r"""
        The scale invariant growth rate :math:`d\mathrm{ln}D/d\mathrm{ln}a` for
        CDM perturbations.

        This is the quantity defined by CLASS as ``index_bg_f`` in the
        background module.
        """
        return self.compute_for_z(z, self.ba.index_bg_f)

cdef class Perturbs:
    """
    A wrapper of the `perturbs module <https://goo.gl/VVhpcU>`_ in CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
    cdef ClassEngine engine
    cdef perturbs * pt
    cdef background * ba

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("perturbs")
        self.pt = &self.engine.pt
        self.ba = &self.engine.ba

    property k_max_for_pk:
        r"""
        The input parameter specifying the maximum ``k`` value to compute
        spectra for; units of :math:`h \mathrm{Mpc}^{-1}`.
        """
        def __get__(self):
            return self.pt.k_max_for_pk/self.ba.h

    property P_z_max:
        r"""
        The input parameter specifying the maximum redshift measured for
        power spectra.
        """
        def __get__(self):
            return self.pt.z_max_pk

    property gauge:
        r"""
        The gauge name as a string, either 'newtonian' or 'synchronous'.
        """
        def __get__(self):
            if self.pt.gauge == newtonian:
              return 'newtonian'
            elif self.pt.gauge == synchronous:
              return 'synchronous'
            else:
              raise ValueError("gauge value not understood")

cdef class Thermo:
    """
    A wrapper of the `thermo module <https://goo.gl/JKGUP6>`_ in CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
    cdef ClassEngine engine
    cdef thermo * th
    cdef background * ba

    def __init__(self, ClassEngine engine):
        self.engine = engine
        self.engine.compute("thermodynamics")
        self.th = &self.engine.th
        self.ba = &self.engine.ba

    property z_drag:
        r"""
        The baryon drag redshift.
        """
        def __get__(self):
            return self.th.z_d

    property rs_drag:
        r"""
        The comoving sound horizon at baryon drag, in :math:`\mathrm{Mpc}/h`.
        """
        def __get__(self):
            return self.th.rs_d * self.ba.h

    property tau_reio:
        r"""
        The reionization optical depth.
        """
        def __get__(self):
            return self.th.tau_reio

    property z_reio:
        r"""
        The reionization redshift.
        """
        def __get__(self):
            return self.th.z_reio

    property z_rec:
        r"""
        The redshift at which the visibility reaches its maximum; equals
        the recombination redshift.
        """
        def __get__(self):
            return self.th.z_rec

    property rs_rec:
        r"""
        The comoving sound horizon at recombination, :math:`z=z_\mathrm{rec}`.
        Units of :math:`\mathrm{Mpc}/h`.
        """
        def __get__(self):
            return self.th.rs_rec * self.ba.h

    property theta_s:
        r"""
        The sound horizon angle at recombination, equal to
        :math:`r_s(z_\mathrm{rec}) / D_a(z_\mathrm{rec})`.
        """
        def __get__(self):
            return self.th.rs_rec / self.th.ra_rec


cdef class Primordial:
    """
    A wrapper of the `primordial module <https://goo.gl/SmxLQz>`_ in CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
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

    def get_pkprim(self, k):
        r"""
        The primoridal spectrum of curvation perturabtion at ``k``, generated by 
        inflation. This is defined as:

        .. math::

            \mathcal{P_R}(k) = A_s \left (\frac{k}{k_0} \right )^{n_s - 1 + 0.5 \ln(k/k_0) (dn_s / d\ln k)}

        See also: equation 2 of `this reference <https://arxiv.org/abs/1303.5076>`_.

        Parameters
        ----------
        k : array_like
          wavenumbers in :math:`h \mathrm{Mpc}^{-1}` units.

        Returns
        -------
        array_like :
          the primordial power
        """
        #generate a new output array of the correct shape by broadcasting input arrays together
        k = np.float64(k) * self.ba.h # convert to 1/Mpc
        out = np.empty(np.broadcast(k).shape, np.float64)

        #generate the iterator over the input and output arrays, does the same thing as
        cdef np.broadcast it = np.broadcast(k,  out)
        cdef int index_md = 0

        while np.PyArray_MultiIter_NOTDONE(it):

                #PyArray_MultiIter_DATA is used to access the pointers the iterator points to
                aval = (<double*>np.PyArray_MultiIter_DATA(it, 0))[0]

                if aval == 0: # forcefully set k == 0 to zero.
                    (<double*>(np.PyArray_MultiIter_DATA(it, 1)))[0] = 0.
                else:
                    if _FAILURE_ == primordial_spectrum_at_k(self.pm, index_md, linear, aval,
                        <double*>(np.PyArray_MultiIter_DATA(it, 1))):
                        raise ClassRuntimeError(self.pm.error_message.decode())

                #PyArray_MultiIter_NEXT is used to advance the iterator
                np.PyArray_MultiIter_NEXT(it)

        # Watch out: no transformation here
        return out

    def get_primordial(self):
        r"""
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.

        The 'output' parameter must be set to something, e.g. 'tCl'.

        Returns
        -------
        array_like :
          structured array containing k-vector and primordial scalar and tensor :math:`P(k)`.
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
    """
    A wrapper of the `spectra module <https://goo.gl/EMti1s>`_ in CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
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
        r"""
        Boolean flag specifying whether the power spectrum is nonlinear.
        """
        def __get__(self):
          return self.nl.method > 0

    property has_pk_matter:
        r"""
        Boolean flag specifying whether matter power spectra have been
        requested as output.
        """
        def __get__(self):
          return self.pt.has_pk_matter

    property P_k_min:
        r"""
        The minimum ``k`` value for which power spectra have been computed in
        :math:`h \mathrm{Mpc}^{-1}`.

        This is computed from the ``ln_k`` array of the Spectra module.
        """
        def __get__(self):
            # factor of 1.001 to avoid bounds errors due to rounding errors
            return 1.001*np.exp(self.sp.ln_k[0])/self.ba.h;

    property P_k_max:
        r"""
        The maximum ``k`` value measured for power spectra in
        :math:`h \mathrm{Mpc}^{-1}`.
        """
        def __get__(self):
            # factor of 0.999 to avoid bounds errors due to rounding errors
            return 0.999*np.exp(self.sp.ln_k[self.sp.ln_k_size-1])/self.ba.h;

    property sigma8:
        r"""
        The amplitude of matter fluctuations at :math:`z=0`.
        """
        def __get__(self):
            return self.sp.sigma8

    property A_s:
        r"""
        The scalar amplitude of the primordial power spectrum at
        :math:`k_\mathrm{pivot}`.
        """
        def __get__(self):
            return self.pm.A_s

    property ln_1e10_A_s:
        r"""
        Return :math:`\log(10^{10}A_s)`.
        """
        def __get__(self):
            return np.log(1e10*self.A_s)

    property n_s:
        r"""
        The tilt of the primordial power spectrum.
        """
        def __get__(self):
            return self.pm.n_s

    property k_pivot:
        r"""
        The primordial power spectrum pivot scale, where the primordial power
        is equal to :math:`A_s`. Units of :math:`h \mathrm{Mpc}^{-1}`.
        """
        def __get__(self):
            return self.pm.k_pivot / self.ba.h

    def sigma8_z(self, z):
        r"""
        Return :math:`\sigma_8(z)`.
        """
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
        r"""
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'dCl' and 'vCl' in the list of
        'output'. The transfer functions can also be computed at higher
        redshift ``z`` provided that 'z_pk' has been set and that ``z`` is
        inside the region spanned by 'z_pk'.

        This function is not vectorized because it returns a vector when
        'ic_size' is greater than 1, and we don't understand 'ic_size'.

        Parameters
        ----------
        z  : float
          redshift (default = 0)
        output_format  : ('class' or 'camb')
          Format transfer functions according to CLASS convention (default)
          or CAMB convention.

        Returns
        -------
        tk : array_like
          array containing transfer functions. ``k`` here is in units of
          :math:`h \mathrm{Mpc}^{-1}`.
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
        r"""
        Gives the pk for a given ``k`` and ``z`` (will be nonlinear if requested
        by the user, linear otherwise).

        .. note::

            There is an additional check to verify if output contains `mPk`,
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
        r"""
        The primary power spectrum result (nonlinear if enabled) on ``k`` and
        ``z`` array.

        Parameters
        ----------
        k : float, array_like
          the wavenumber in units of :math:`h \mathrm{Mpc}^{-1}`
        z : float, array_like
          the redshift values

        Returns
        -------
        array like :
            the power spectrum in units of :math:`(\mathrm{Mpc}/h)^3`
        """
        return self._get_pk(k, z, 0)

    def get_pklin(self, k, z):
        r"""
        Linear power spectrum result (linear even if nonlinear is enabled)
        on ``k`` and ``z`` array.

        Parameters
        ----------
        k : float, array_like
          the wavenumber in units of :math:`h \mathrm{Mpc}^{-1}`
        z : float, array_like
          the redshift values

        Returns
        -------
        array like :
            the power spectrum in units of :math:`(\mathrm{Mpc}/h)^3`
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
