from .binding import ClassEngine, Background, Spectra, Perturbs, Primordial
import numpy
from six import string_types
import os


aliases = {'Omega_b': 'Omega0_b', 'Omega_cdm':'Omega0_cdm'}

conflicts = {'h': ['H0', '100*theta_s'],
             'T_cmb': ['Omega_g', 'omega_g'],
             'Omega_b': ['omega_b'],
             'N_ur': ['Omega_ur', 'omega_ur'],
             'Omega_cdm': ['omega_cdm'],
             'm_ncdm': ['Omega_ncdm', 'omega_ncdm'],
             'P_k_max': ['P_k_max_h/Mpc', 'P_k_max_1/Mpc'],
             'P_z_max': ['z_max_pk'],
             'sigma8': ['A_s', 'ln10^{10}A_s'],
             'nonlinear' : ['non linear']
            }

class Cosmology(object):
    """
    A cosmology calculator based on the CLASS binding in classylss.

    It is a collection of all method provided by the class interfaces.

    The individual interfaces can be accessed too, such that
    `c.Spectra.get_transfer` and `c.get_transfer` are identical.

    Notes
    -----
    * The default configuration assumes a flat cosmology, :math:`\Omega_{0,k}=0`.
      Pass ``Omega0_k`` in the ``extra`` keyword dictionary to change this value.
    * By default, a cosmological constant is assumed, with its density value
      inferred by the curvature condition.
    * Non-cosmological constant dark energy can be used by specifying the
      ``w0_fld``, ``wa_fld``, and/or ``Omega0_fld`` values.

    Parameters
    ----------
    h : float
        the dimensionaless Hubble parameter
    Tcmb0 : float
        the temperature of the CMB in Kelvins
    Omega0_b : float
        the current baryon density parameter, :math:`\Omega_{b,0}`
    Omega0_cdm : float
        the current cold dark matter density parameter, :math:`\Omega_{cdm,0}`
    N_ur : float
        the number of ultra-relativistic (massless neutrino) species; the
        default number is inferred based on the number of massive neutrinos
        via the following logic: if you have respectively 1,2,3 massive
        neutrinos and use the default ``T_ncdm`` value (0.71611 K), designed
        to give m/omega of 93.14 eV, and you wish to have ``N_eff=3.046`` in
        the early universe, then ``N_ur`` is set to 2.0328, 1.0196, 0.00641,
        respectively.
    m_ncdm : list
        the masses (in eV) for all massive neutrino species; an empty list
        should  be passed for no massive neutrinso. The default is a single
        massive neutrino with mass of 0.06 eV
    P_k_max : float
        the maximum ``k`` value to compute power spectrum results to, in units
        of :math:`h/Mpc`
    P_z_max : float
        the maximum redshift to compute power spectrum results to
    gauge : str,
        either synchronous or newtonian
    sigma8 : float
        the desired current amplitude of matter fluctuations; the scalar
        amplitude parameter ``A_s`` will be automatically adjusted to
        achieve the desired ``sigma8`` value
    n_s : float
        the tilt of the primordial power spectrum
    nonlinear : bool
        whether to compute nonlinear power spectrum results via HaloFit
    verbose : bool
        whether to turn on the default CLASS logging for all submodules
    extra : dict
        dictionary of extra parameters to pass to CLASS; users should be wary
        of configuration options that may conflict with the base set
        of parameters
    """
    # delegate resolve order -- a pun at mro; which in
    # this case introduces the meta class bloat and doesn't solve
    # the issue. We want delayed initialization of interfaces
    # or so-called 'mixins'.
    # easier to just use delegates with a customized getattr.
    # this doesn't work well with automated documentation tools though,
    # unfortunately.

    dro = [Spectra, Perturbs, Primordial, Background, ClassEngine]
    dro_dict = dict([(n.__name__, n) for n in dro])

    def __init__(self,
            h=0.67556,
            T_cmb=2.7255,
            Omega_b=0.022032/0.67556**2,
            Omega_cdm=0.12038/0.67556**2,
            N_ur=None,
            m_ncdm=[0.06],
            P_k_max=10.,
            P_z_max=100.,
            gauge='synchronous',
            sigma8=0.82,
            n_s=0.9667,
            nonlinear=False,
            verbose=False,
            extra={}, # additional arguments to pass to CLASS
            **kwargs
        ):
        # quickly copy over all arguments --
        # at this point locals only contains the arguments.
        args = dict(locals())

        # store the extra CLASS params
        extra = args.pop('extra')
        extra.update(args.pop('kwargs'))

        # remove some non-CLASS variables
        args.pop('self')

        # use set state to de-serialize the object.
        self.__setstate__((args,extra))

    def __dir__(self):
        """ a list of all members from all delegate classes """
        r = []
        # first allow tab completion of delegate names; to help resolve conflicts
        r.extend([n.__name__ for n in self.dro])
        # then allow tab completion of all delegate methods
        for i in reversed(self.dro):
            r.extend(dir(i))
        return sorted(list(set(r)))

    def __getattr__(self, name):
        """ Will find the proper delegate, initialize it, and run the method """
        # getting a delegate explicitly, e.g. c.Background
        if name in self.dro_dict:
            iface = self.dro_dict[name]
            if iface not in self.delegates:
                self.delegates[iface] = iface(self.engine)
            return self.delegates[iface]

        # resolving a name from the delegates : c.Om0 => c.Background.Om0
        for iface in self.dro:
            if hasattr(iface, name):
                if iface not in self.delegates:
                    self.delegates[iface] = iface(self.engine)
                d = self.delegates[iface]
                return getattr(d, name)
        else:
            raise AttributeError("Attribute `%s` not found in any of the delegate objects" % name)

    def __getstate__(self):
        return (self.args, self.extra)

    @property
    def sigma8(self):
        """
        The present day value of ``sigma_r(r=8 Mpc/h)``, used to normalize
        the power spectrum, which is proportional to the square of this value.

        The power spectrum can re-normalized by setting a different
        value for this parameter
        """
        return self.Spectra.sigma8

    @sigma8.setter
    def sigma8(self, value):
        """
        Set the sigma8 value and normalize the power spectrum to the new value
        """
        if not numpy.isclose(self.sigma8, value):
            extra = self.extra.copy()
            for key in conflicts['sigma8']:
                if key in extra: extra.pop(key)
            A_s = self.A_s * (value/self.sigma8)**2
            extra['A_s'] = A_s
            self.__setstate__((self.args,extra))

    @classmethod
    def from_file(cls, filename):
        """
        Initialize a :class:`Cosmology` object from the CLASS parameter file

        Parameters
        ----------
        filename : str
            the name of the parameter file to read
        """
        # make sure it is a valid file
        if not os.path.exists(filename):
            raise ValueError("no such file: %s" %filename)

        # extract dictionary of parameters from the file
        fc = open(filename, 'r').read()
        pars = read_CLASS_ini(filename)

        # initialize the engine as the backup delegate.
        toret = object.__new__(cls)
        toret.engine = ClassEngine(pars)
        toret.delegates = {ClassEngine: toret.engine}

        # reconstruct the correct __init__ params
        args = {}
        for name in list(pars.keys()):
            val = pars.pop(name)

            if name in conflicts:
                alias = aliases.get(name, name)
                args[name] = getattr(toret, alias)
            else:
                for c in conflicts:
                    if c == 'sigma8': continue
                    if name in conflicts[c]:
                        alias = aliases.get(c, c)
                        args[c] = getattr(toret, alias)

        # set the gauge
        args['gauge'] = 'newtonian' if toret.gauge == 0 else 'synchronous'

        # set sigma8 norm
        if not any(par in pars for par in conflicts['sigma8']):
            pars['A_s'] = toret.A_s

        toret.args = args
        toret.extra = pars

        return toret

    def __setstate__(self, state):

        # remember for serialization
        self.args, self.extra = state

        # remove sigma8 for norm later
        desired_sigma8 = self.args.pop('sigma8', None)

        # verify and set defaults
        pars = verify_parameters(self.args, self.extra)

        # initialize the engine as the backup delegate.
        self.engine = ClassEngine(pars)
        self.delegates = {ClassEngine: self.engine}

        # set sigma8 by re-scaling A_s, if A_s was not specified
        if desired_sigma8 is not None:
            if not any(key in self.extra for key in conflicts['sigma8']):
                self.sigma8 = desired_sigma8

    def clone(self, **kwargs):
        """
        Create a new cosmology based on modification of self
        """
        desired_sigma8 = kwargs.pop('sigma8', None)
        if desired_sigma8 is not None:
            kwargs['A_s'] = self.A_s * (desired_sigma8/self.sigma8)**2

        new = Cosmology(**kwargs)
        args = self.args.copy()
        args.update(new.args)

        args['extra'] = self.extra.copy()
        args['extra'].update(new.extra)

        return Cosmology(**args)

    @classmethod
    def from_astropy(kls, cosmo, extra={}):
        args = astropy_to_dict(cosmo)
        args.update(extra)
        return Cosmology(**args)

def astropy_to_dict(cosmo):

    from astropy import cosmology, units

    pars = {'extra':{}}
    pars['h'] = cosmo.h
    pars['T_cmb'] = cosmo.Tcmb0
    if cosmo.Ob0 is not None:
        pars['Omega0_b'] = cosmo.Ob0
    else:
        raise ValueError("please specify a value 'Ob0' ")
    pars['Omega0_cdm'] = cosmo.Om0 - cosmo.Ob0 # should be okay for now

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

        pars['m_ncdm'] = [k.value for k in sorted(m_nu[m_nu > 0.], reverse=True)]
    else:
        pars['m_ncdm'] = []
        pars['N_ur'] = cosmo.Neff

    # specify the curvature
    pars['extra']['Omega_k'] = cosmo.Ok0

    # handle dark energy
    if isinstance(cosmo, cosmology.LambdaCDM):
        pass
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
        raise ValueError(msg)

    return pars

def read_CLASS_ini(filename):
    """
    Read a CLASS ``.ini`` file, returning a dictionary of parameters

    Parameters
    ----------
    filename : str
        the name of the file to read
    """
    pars = {}

    with open(filename, 'r') as ff:

        # loop over lines
        for lineno, line in enumerate(ff):

            # skip any commented lines with #
            if '#' in line: line = line[line.index('#'):]

            # must have an equals sign to be valid
            if "=" not in line: continue

            # extract key and value pairs
            fields = line.split("=")
            assert len(fields) == 2, "error reading line number %d" %lineno
            pars[fields[0].strip()] = fields[1].strip()

    return pars

def verify_parameters(args, extra):
    """
    Verify the input parameters to a :class:`Cosmology` object and
    set various default values
    """
    # check for conflicts
    for par in conflicts:
        for p in conflicts[par]:
            if p in extra and par is not 'sigma8':
                raise ValueError("input parameter conflict; use '%s', not '%s'" %(par, p))

    pars = {}
    pars.update(args)
    pars.update(extra)

    # set some default parameters
    pars.setdefault('output', "vTk dTk mPk")
    pars.setdefault('extra metric transfer functions', 'y')

    # no massive neutrinos
    if pars.get('m_ncdm', None) is None:
        pars['m_ncdm'] = []

    # a single massive neutrino
    if numpy.isscalar(pars['m_ncdm']):
        pars['m_ncdm'] = [pars['m_ncdm']]

    # needs to be a list
    if not isinstance(pars['m_ncdm'], list):
        raise TypeError("``m_ncdm`` should be a list of mass values in eV")

    # check gauge
    if pars.get('gauge', 'synchronous') not in ['synchronous', 'newtonian']:
        raise ValueError("'gauge' should be 'synchronous' or 'newtonian'")

    for m in pars['m_ncdm']:
        if m == 0:
            raise ValueError("A zero mass is specified in the massive neutrino list.")

    # remove None's -- use CLASS default
    for key in list(pars.keys()):
        if pars[key] is None: pars.pop(key)

    # turn on verbosity
    verbose = pars.pop('verbose', False)
    if verbose:
        for par in ['input', 'background', 'thermodynamics', 'perturbations',
                    'transfer', 'primordial', 'spectra', 'nonlinear', 'lensing']:
            name = par + '_verbose'
            if name not in pars: pars[name] = 1

    # maximum k value
    if 'P_k_max_h/Mpc' not in pars:
        pars['P_k_max_h/Mpc'] = pars.pop('P_k_max', 10.)

    # maximum redshift
    if 'z_max_pk' not in pars:
        pars['z_max_pk'] = pars.pop('P_z_max', 100.)

    # nonlinear power?
    if 'non linear' not in pars:
        if pars.pop('nonlinear', False):
            pars['non linear'] = 'halofit'

    # from CLASS notes:
    # one more remark: if you have respectively 1,2,3 massive neutrinos,
    # if you stick to the default value pm equal to 0.71611, designed to give m/omega of
    # 93.14 eV, and if you want to use N_ur to get N_eff equal to 3.046 in the early universe,
    # then you should pass here respectively 2.0328,1.0196,0.00641
    N_ur_table = [3.046, 2.0328, 1.0196, 0.00641]
    if 'N_ur' not in pars:
        pars['N_ur'] = N_ur_table[len(pars['m_ncdm'])]
    pars['N_ncdm'] = len(pars['m_ncdm'])

    return pars
