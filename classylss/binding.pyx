#cython: embedsignature=True
cimport cython
import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset, strncpy, strdup

from classylss import get_data_files

_DATA_FILES = get_data_files()

from cclassy cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass

cdef void _build_file_content(pars, file_content * fc):
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
    for kk in pars:

        dumcp = kk.encode()
        strncpy(fc.name[i], dumcp[:sizeof(FileArg)-1], sizeof(FileArg))
        dumcp = str(pars[kk]).encode()
        strncpy(fc.value[i], dumcp[:sizeof(FileArg)-1], sizeof(FileArg))
        fc.read[i] = _FALSE_

        i += 1

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


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass

ctypedef struct ready_flags:
    int ba
    int th
    int pt
    int pm
    int nl
    int tr
    int sp
    int op
    int le


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

    def __cinit__(self, *args, **kwargs):
        memset(&self.ready, 0, sizeof(self.ready))

    def __init__(self, object pars={}, level=['lensing']):
        cdef char* dumc

        cdef file_content fc

        tasks = _build_task_dependency(level)

        _build_file_content(pars, &fc)

        try:
            self.compute(&fc, tasks)
        finally:
            # free fc regardless of compute exceptions.
            parser_free(&fc) 
            pass

    def __dealloc__(self):
        if self.ready.le: lensing_free(&self.le)
        if self.ready.sp: spectra_free(&self.sp)
        if self.ready.tr: transfer_free(&self.tr)
        if self.ready.nl: nonlinear_free(&self.nl)
        if self.ready.pm: primordial_free(&self.pm)
        if self.ready.pt: perturb_free(&self.pt)
        if self.ready.th: thermodynamics_free(&self.th)
        if self.ready.ba: background_free(&self.ba)

    cdef compute(self, file_content * fc, object tasks):
        """
        Main function, execute all the _init methods for all desired modules.
        This is called in MontePython, and this ensures that the Class instance
        of this class contains all the relevant quantities. Then, one can deduce
        Pk, Cl, etc...

        Parameters
        ----------
        level : list
                list of the last module desired. The internal function
                _check_task_dependency will then add to this list all the
                necessary modules to compute in order to initialize this last
                one. The default last module is "lensing".

        .. warning::

            level default value should be left as an array (it was creating
            problem when casting as a set later on, in _check_task_dependency)

        """
        cdef ErrorMsg errmsg

        level = tasks

        # --------------------------------------------------------------------
        # Check the presence for all CLASS modules in the list 'level'. If a
        # module is found in level, executure its "_init" method.
        # --------------------------------------------------------------------
        # The input module should raise a CosmoSevereError, because
        # non-understood parameters asked to the wrapper is a problematic
        # situation.
        if "input" in level:
            if input_init(fc, &self.pr, &self.ba, &self.th,
                          &self.pt, &self.tr, &self.pm, &self.sp,
                          &self.nl, &self.le, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            for i in range(fc.size):
                if fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(fc.name[i].decode())
            if problem_flag:
                raise CosmoSevereError(
                    "Class did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a CosmoComputationError
        # with the error message from the faulty module of CLASS.
        if "background" in level:
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                raise CosmoComputationError(self.ba.error_message)
            self.ready.ba = True

        if "thermodynamics" in level:
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                raise CosmoComputationError(self.th.error_message)
            self.ready.th = True

        if "perturb" in level:
            if perturb_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                raise CosmoComputationError(self.pt.error_message)
            self.ready.pt = True

        if "primordial" in level:
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                raise CosmoComputationError(self.pm.error_message)
            self.ready.pm = True

        if "nonlinear" in level:
            if nonlinear_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.nl) == _FAILURE_:
                raise CosmoComputationError(self.nl.error_message)
            self.ready.nl = True

        if "transfer" in level:
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.nl), &(self.tr)) == _FAILURE_:
                raise CosmoComputationError(self.tr.error_message)
            self.ready.tr = True

        if "spectra" in level:
            if spectra_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.nl), &(self.tr),
                            &(self.sp)) == _FAILURE_:
                raise CosmoComputationError(self.sp.error_message)
            self.ready.sp = True

        if "lensing" in level:
            if lensing_init(&(self.pr), &(self.pt), &(self.sp),
                            &(self.nl), &(self.le)) == _FAILURE_:
                raise CosmoComputationError(self.le.error_message)
            self.ready.le = True

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return


    @classmethod
    def from_astropy(cls, cosmo, extra={}):
        """
        Convert an astropy cosmology to a ``ClassParams`` instance
        """
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
            raise TypeError(msg)

        # add any extra arguments
        if len(extra):
            pars.update(extra)

        _pars = {
            "output": "vTk dTk mPk",
            "P_k_max_h/Mpc":  20.,
            "z_max_pk": 100.0,
            }

        for k in _pars:
            pars.setdefault(k, _pars[k])

        self = cls(pars)
        return self

cdef class Background:
    cdef ClassEngine engine
    cdef background * ba
    cdef readonly dict data

    def __init__(self, ClassEngine engine):
        self.engine = engine
        assert engine.ready.ba
        self.ba = &self.engine.ba
        self.data = self._get_background()

    def Hubble(self, z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')

        if background_tau_of_z(self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,&pvecback[0])==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        return H

    property Omega_m:
        def __get__(self):
            return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm

    property Omega_b:
        def __get__(self):
            return self.ba.Omega0_b

    property omega_b:
        def __get__(self):
            return self.ba.Omega0_b * self.ba.h * self.ba.h

    property Omega_nu:
        def __get__(self):
            return self.ba.Omega0_ncdm_tot

    property Neff:
        def __get__(self):
            return self.ba.Neff

    property age:
        def __get__(self):
            return self.ba.age

    property h:
        def __get__(self):
            return self.ba.h

    property T_cmb:
        """
        Return the CMB temperature
        """
        def __get__(self):
            return self.ba.T_cmb

    property Omega0_m:
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        def __get__(self):
            return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm


    def z_of_r (self, z_array):
        cdef double tau=0.0
        cdef int last_index=0 #junk
        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')
        r = np.zeros(len(z_array),'float64')
        dzdr = np.zeros(len(z_array),'float64')

        i = 0
        for redshift in z_array:
            if background_tau_of_z(self.ba,redshift,&tau)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            if background_at_tau(self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,&pvecback[0])==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # store r
            r[i] = pvecback[self.ba.index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[self.ba.index_bg_H]

            i += 1

        return r[:],dzdr[:]

    def luminosity_distance(self, z):
        """
        luminosity_distance(z)
        """
        cdef double tau=0.0
        cdef int last_index = 0  # junk
        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')

        if background_tau_of_z(self.ba, z, &tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(self.ba, tau, self.ba.long_info,
                self.ba.inter_normal, &last_index, &pvecback[0])==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)
        lum_distance = pvecback[self.ba.index_bg_lum_distance]
        return lum_distance

    def angular_distance(self, z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')

        if background_tau_of_z(self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,&pvecback[0])==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        return D_A

    def scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)

        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double [::1] pvecback = np.zeros(self.ba.bg_size, dtype='f8')

        if background_tau_of_z(self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index, &pvecback[0])==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        return D

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
        assert engine.ready.ba
        assert engine.ready.sp
        assert engine.ready.pt
        assert engine.ready.pm
        assert engine.ready.nl

        self.ba = &self.engine.ba
        self.nl = &self.engine.nl
        self.sp = &self.engine.sp
        self.pt = &self.engine.pt
        self.pm = &self.engine.pm

    def get_transfer(self, z=0., output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'dCl' and 'vCl' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that z is inside the region spanned by 'z_pk'.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : dictionary containing transfer functions.
        """

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            raise RuntimeError("Perturbation is not computed")

        cdef FileName ic_suffix
        cdef file_format_outf
        cdef char ic_info[1024]
        cdef char titles[_MAXTITLESTRINGLENGTH_]

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        memset(titles, 0, _MAXTITLESTRINGLENGTH_)

        if spectra_output_tk_titles(self.ba, self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.op.error_message)

        tmp = (<bytes>titles).decode()
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)

        index_md = 0
        ic_num = self.sp.ic_size[index_md]

        dtype = [(name.split()[0], 'f8') for name in names]

        cdef np.ndarray data = np.zeros((ic_num, self.sp.ln_k_size), dtype=dtype)

        if spectra_output_tk_data(self.ba, self.pt, self.sp, outf, <double> z, number_of_titles, <double*> data.data)==_FAILURE_:
            raise CosmoSevereError(self.sp.error_message)

        ic_keys = []
        if ic_num > 1:
            for index_ic in range(ic_num):
                if spectra_firstline_and_ic_suffix(self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                    raise CosmoSevereError(self.op.error_message)

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
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        if not lin:
            if (self.nl.method == 0):
                 if spectra_pk_at_k_and_z(self.ba,self.pm,self.sp,k,z,pk,pk_ic)==_FAILURE_:
                     raise CosmoSevereError(self.sp.error_message)
            else:
                 if spectra_pk_nl_at_k_and_z(self.ba,self.pm,self.sp,k,z,pk) ==_FAILURE_:
                        raise CosmoSevereError(self.sp.error_message)

        else:
            if spectra_pk_at_k_and_z(self.ba,self.pm,self.sp,k,z,pk,pk_ic)==_FAILURE_:
                raise CosmoSevereError(self.sp.error_message)

        return 0

    def get_pk(self, k, z):
        """ Fast function to get the power spectrum on a k and z array """
        from numpy.lib.stride_tricks import broadcast_arrays
        k1, z1 = np.float64(k), np.float64(z)
        return self._get_pk(k1, z1, 0)

    def get_pklin(self, k, z):
        """ Fast function to get the power spectrum on a k and z array """
        from numpy.lib.stride_tricks import broadcast_arrays
        k1, z1 = np.float64(k), np.float64(z)
        return self._get_pk(k1, z1, 1)

    cdef _get_pk(self, k, z, int linear):

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )

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

        return out

