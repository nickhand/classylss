from .binding import Background
import numpy as np

class AstropyCompat(object):
    """ A cosmology delegate object that provides astropy compatibility"""
    def __init__(self, engine):
        self.engine = engine
        self.bg = Background(engine)

    @property
    def Om0(self): return self.bg.Omega0_m
    def Om(self, z): return self.bg.Omega_m(z)

    @property
    def Ob0(self): return self.bg.Omega0_b
    def Ob(self, z): return self.bg.Omega_b(z)

    @property
    def Ogamma0(self): return self.bg.Omega0_g
    def Ogamma(self, z): return self.bg.Omega_g(z)

    @property
    def Odm0(self): return self.bg.Omega0_cdm
    def Odm(self, z): return self.bg.Omega_cdm(z)

    @property
    def Ok0(self): return self.bg.Omega0_k
    def Ok(self, z): return self.bg.Omega_k(z)

    @property
    def Ode0(self): return self.bg.Omega0_lambda + self.bg.Omega0_fld
    def Ode(self, z): return self.bg.Omega_lambda(z) + self.bg.Omega_fld(z)

    @property
    def Onu0(self):
        return self.bg.Omega0_ncdm_tot + self.bg.Omega0_ur
    def Onu(self, z):
        return self.bg.Omega_ncdm(z) + self.bg.Omega_ur(z)

    @property
    def Tcmb0(self): return self.bg.T0_cmb
    def Tcmb(self, z): return self.bg.T0_cmb * (1 + z)

    @property
    def Tnu0(self): return self.bg.T0_ncdm
    def Tnu(self, z): return self.bg.T0_ncdm * (1 + z)

    @property
    def w0(self): return self.bg.w0_fld
    @property
    def wa(self): return self.bg.wa_fld

    @property
    def has_massive_nu(self):
        return self.bg.N_ncdm > 0

    def nu_relative_density(self, z):
        return self.Onu(z) / self.Ogamma(z)
