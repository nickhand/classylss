from .binding import Background
import numpy as np

class AstropyCompat(object):
    """
    A cosmology wrapper object that provides compatibility with the syntax
    from :mod:`astropy`.

    .. warning ::
        Some parameter definitions vary between :mod:`astropy` and CLASS.

    Parameters
    ----------
    engine : ClassEngine
      the CLASS engine object
    """
    def __init__(self, engine):
        self.engine = engine
        self.bg = Background(engine)

    @property
    def Om0(self):
        """Returns :attr:`~classylss.binding.Background.Omega0_m`."""
        return self.bg.Omega0_m
    def Om(self, z):
        """Returns :func:`~classylss.binding.Background.Omega_m`."""
        return self.bg.Omega_m(z)

    @property
    def Ob0(self):
        """Returns :attr:`~classylss.binding.Background.Omega0_b`."""
        return self.bg.Omega0_b
    def Ob(self, z):
        """Returns :func:`~classylss.binding.Background.Omega_b`."""
        return self.bg.Omega_b(z)

    @property
    def Ogamma0(self):
        """Returns :attr:`~classylss.binding.Background.Omega0_g`."""
        return self.bg.Omega0_g
    def Ogamma(self, z):
        """Returns :func:`~classylss.binding.Background.Omega_g`."""
        return self.bg.Omega_g(z)

    @property
    def Odm0(self):
        """Returns :attr:`~classylss.binding.Background.Omega0_cdm`."""
        return self.bg.Omega0_cdm
    def Odm(self, z):
        """Returns :func:`~classylss.binding.Background.Omega_cdm`."""
        return self.bg.Omega_cdm(z)

    @property
    def Ok0(self):
        """Returns :attr:`~classylss.binding.Background.Omega0_k`."""
        return self.bg.Omega0_k
    def Ok(self, z):
        """Returns :func:`~classylss.binding.Background.Omega_k`."""
        return self.bg.Omega_k(z)

    @property
    def Ode0(self):
        """
        Returns the sum of :attr:`~classylss.binding.Background.Omega0_lambda`
        and :attr:`~classylss.binding.Background.Omega0_fld`.
        """
        return self.bg.Omega0_lambda + self.bg.Omega0_fld
    def Ode(self, z):
        """
        Returns the sum of :func:`~classylss.binding.Background.Omega_lambda`
        and :func:`~classylss.binding.Background.Omega_fld`.
        """
        return self.bg.Omega_lambda(z) + self.bg.Omega_fld(z)

    @property
    def Onu0(self):
        """
        Returns the sum of :attr:`~classylss.binding.Background.Omega0_ncdm_tot`
        and :attr:`~classylss.binding.Background.Omega0_ur`.
        """
        return self.bg.Omega0_ncdm_tot + self.bg.Omega0_ur
    def Onu(self, z):
        """
        Returns the sum of :func:`~classylss.binding.Background.Omega_ncdm`
        and :func:`~classylss.binding.Background.Omega_ur`.
        """
        return self.bg.Omega_ncdm(z) + self.bg.Omega_ur(z)

    @property
    def Tcmb0(self):
        """Returns :attr:`~classylss.binding.Background.T0_cmb`."""
        return self.bg.T0_cmb
    def Tcmb(self, z):
        """Returns :math:`(1+z)` :attr:`~classylss.binding.Background.T0_cmb`."""
        return self.bg.T0_cmb * (1 + z)

    @property
    def Tnu0(self):
        """Returns :attr:`~classylss.binding.Background.T0_ncdm`."""
        return self.bg.T0_ncdm
    def Tnu(self, z):
        """Returns :math:`(1+z)` :attr:`~classylss.binding.Background.T0_ncdm`."""
        return self.bg.T0_ncdm * (1 + z)

    @property
    def w0(self):
        """Returns :attr:`~classylss.binding.Background.w0_fld`."""
        return self.bg.w0_fld
    @property
    def wa(self):
        """Returns :attr:`~classylss.binding.Background.wa_fld`."""
        return self.bg.wa_fld

    @property
    def has_massive_nu(self):
        """
        Returns ``True`` if :attr:`~classylss.binding.Background.N_ncdm`
        is greater than zero.
        """
        return self.bg.N_ncdm > 0

    def nu_relative_density(self, z):
        """
        Returns :func:`Onu`  / :func:`Ogamma`.
        """
        return self.Onu(z) / self.Ogamma(z)
