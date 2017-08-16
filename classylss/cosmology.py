from .binding import ClassEngine, Background, Spectra, Perturbs, Primordial

class Cosmology(object):
    """
        A cosmology calculator based on the CLASS binding in classylss

        It is a collection of all method provided by the class interfaces.

        The individual interfaces can be accessed too, such that
        `c.Spectra.get_transfer` and `c.get_transfer` are identical.

        Parameters are not well documented, and may grow / shrink over time.

        'h' : 'Hubble Parameter',
        'T_cmb' : 'photon ',
        'Omega_b' : 'Baryon',
        'Omega_cdm' : 'CDM',
        'N_ur' : 'Number of Ultra relativistic species',
        'N_ncdm' : 'Number of non-relativistic species',
        'm_ncdm' : neutrino masses
        'w0_fld' : 'w0',
        'wa_fld'  : 'wa',
        'Omega_Lambda' : 'Lambda, 0 if w0 is given',
         gauge : synchronous or newtonian

    """

    # delegate resolve order -- a pun at mro; which in
    # this case introduces the meta class bloat and doesn't solve
    # the issue. We want delayed initialization of interfaces
    # or so-called 'mixins'.
    # easier to just use delegates with a customized getattr.
    # this doesn't work well with automated documentation tools though,
    # unfortunately.

    dro = [Spectra, Perturbs, Primordial, Background, ClassEngine]

    def __init__(self,
            h=0.67556,
            T_cmb=2.7255,
            Omega_b=0.022032/0.67556**2,
            Omega_cdm=0.12038/0.67556**2,
            Omega_Lambda=None,
            N_ur=None,
            m_ncdm=[0.06],
            P_k_max=10.,
            P_z_max=100.,
            gauge='synchronous',
            **kwargs # additional arguments from CLASS
        ):
        # quickly copy over all arguments -- 
        # at this point locals only contains the arguments.

        args = dict(locals())
        args.pop('self')
        kwargs = args.pop('kwargs')
        args.update(kwargs)

        # use set state to de-serialize the object.
        self.__setstate__((args,))

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
        if name in self.dro:
            iface = name
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
            raise AttributeError("Attribute `%s` not found in any of the delegate objects")

    def __getstate__(self):
        return (self.args, )

    def __setstate__(self, state):
        args = state[0]

        # remember args for serialization
        self.args = args

        pars = {
            "output": "vTk dTk mPk",
            "extra metric transfer functions": "y",
            }
        pars.update(args)

        # remove None's -- use CLASS default
        for key in list(pars.keys()):
            if pars[key] is None: pars.pop(key)

        assert args['gauge'] in ['synchronous', 'newtonian']

        if 'P_k_max_h/Mpc' not in pars:
            pars['P_k_max_h/Mpc'] = pars.pop('P_k_max')

        if 'z_max_pk' not in pars:
            pars['z_max_pk'] = pars.pop('P_z_max')

        for m in args['m_ncdm']:
            if m == 0:
                raise ValueError("A zero mass is specified in the massive neutrino list.")

        # from CLASS notes:
        # one more remark: if you have respectively 1,2,3 massive neutrinos,
        # if you stick to the default value pm equal to 0.71611, designed to give m/omega of
        # 93.14 eV, and if you want to use N_ur to get N_eff equal to 3.046 in the early universe,
        # then you should pass here respectively 2.0328,1.0196,0.00641

        N_ur_table = [3.046, 2.0328, 1.0196, 0.00641]

        if 'N_ur' not in pars:
            pars['N_ur'] = N_ur_table[len(args['m_ncdm'])]

        pars['N_ncdm'] = len(args['m_ncdm'])

        # initialize the engine as the backup delegate.
        self.engine = ClassEngine(pars)
        self.delegates = {
            ClassEngine: self.engine,
        }

    def clone(self, **kwargs):
        """ create a new cosmology based on modification of self """
        args = {}
        args.update(self.args)
        args.update(self.kwargs)
        return Cosmology(**args)

    @classmethod
    def from_astropy(kls, cosmo, extra={}):
        args = astropy_to_dict(cosmo)
        args.update(extra)
        return Cosmology(**args)

def astropy_to_dict(cosmo):

    from astropy import cosmology, units

    pars = {}
    pars['h'] = cosmo.h
    pars['T_cmb'] = cosmo.Tcmb0
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

        pars['m_ncdm'] = [k.value for k in sorted(m_nu[m_nu > 0.], reverse=True)]
    else:
        pars['m_ncdm'] = []
        pars['N_ur'] = cosmo.Neff

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

    return pars

