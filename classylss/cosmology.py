from .binding import ClassEngine, Background, Spectra, Perturbs, Primordial

class Cosmology(object):
    """
        'h' : 'Hubble Parameter',
        'Omega_g' : 'photon ',
        'Omega_b' : 'Baryon',
        'Omega_cdm' : 'CDM',
        'N_ur' : 'Number of Ultra relativistic species',
        'N_ncdm' : 'Number of non-relativistic species',
        'm_ncdm' : 'comma sperated string of neutrino mass, length of N_ncdm',
        'w0_fld' : 'w0',
        'wa_fld'  : 'wa',
        'Omega_Lambda' : 'Lambda, 0 if w0 is given',
    """
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
            **kwargs
        ):
        # quickly copy over all arguments -- 
        # at this point locals only contains the arguments.

        args = dict(locals())
        args.pop('self')
        kwargs = args.pop('kwargs')

        args.update(kwargs)

        # remove None's -- use CLASS default
        for key in list(args.keys()):
            if args[key] is None: args.pop(key)

        if 'P_k_max_h/Mpc' not in args:
            args['P_k_max_h/Mpc'] = args.pop('P_k_max')

        if 'z_max_pk' not in args:
            args['z_max_pk'] = args.pop('P_z_max')

        pars = {
            "output": "vTk dTk mPk",
            "extra metric transfer functions": "y",
            }
        pars.update(args)

        for m in m_ncdm:
            if m == 0:
                raise ValueError("A zero mass is specified in the massive neutrino list.")

        # from CLASS notes:
        # one more remark: if you have respectively 1,2,3 massive neutrinos,
        # if you stick to the default value pm equal to 0.71611, designed to give m/omega of
        # 93.14 eV, and if you want to use N_ur to get N_eff equal to 3.046 in the early universe,
        # then you should pass here respectively 2.0328,1.0196,0.00641

        N_ur_table = [3.046, 2.0328, 1.0196, 0.00641]

        if 'N_ur' not in pars:
            pars['N_ur'] = N_ur_table[len(m_ncdm)]

        pars['N_ncdm'] = len(m_ncdm)
        self.engine = ClassEngine(pars)

        # FIXME: remove this
        print(self.engine.parameter_file)

        self.delegates = {
            ClassEngine: self.engine,
        }

        # delegate resolve order
        self.dro = [Spectra, Perturbs, Primordial, Background, ClassEngine]

    def __dir__(self):
        r = []
        for i in reversed(self.dro):
            r.extend(dir(i))
        return sorted(list(set(r)))

    def _resolve(self, name):
        for iface in self.dro:
            if hasattr(iface, name):
                if iface not in self.delegates:
                    self.delegates[iface] = iface(self.engine)
                    print(self.delegates[iface])
                return self.delegates[iface]

    def __getattr__(self, name):
        d = self._resolve(name)
        return getattr(d, name)


