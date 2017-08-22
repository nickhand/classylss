import numpy
from . import transfers

class LinearPower(object):
    """
    An object to compute the linear power spectrum and related quantities,
    using a transfer function from the CLASS code or the analytic
    Eisenstein & Hu approximation.

    Parameters
    ----------
    cosmo : :class:`Cosmology`
        the cosmology instance
    z : float
        the redshift of the power spectrum
    transfer : str, optional
        string specifying the transfer function to use; one of
        'CLASS', 'EisensteinHu', 'NoWiggleEisensteinHu'
    """
    def __init__(self, cosmo, z, transfer='CLASS'):

        self.cosmo = cosmo

        # set sigma8 to the cosmology value
        self._sigma8 = self._sigma8_0 = self.cosmo.sigma8

        # setup the transfers
        if transfer not in transfers.available:
            raise ValueError("'transfer' should be one of %s" %str(transfers.available))
        self.transfer = transfer
        self._transfer = getattr(transfers, transfer)(cosmo, z)
        self._fallback = transfers.EisensteinHu(cosmo, z) # fallback to analytic when out of range

        # normalize to proper sigma8
        self._norm = 1.
        self.z = 0; s8_0 = self.sigma_r(8.) # sigma_r(z=0, r=8)
        self.z = z; s8_z = self.sigma_r(8.) # sigma_r(z=z, r=8)
        self._norm =  (self._sigma8/s8_0)**2 * (s8_z / s8_0)**2

    @property
    def z(self):
        """
        The redshift of the power spectrum
        """
        return self._z

    @z.setter
    def z(self, value):
        self._z = value
        self._transfer.z = value
        self._fallback.z = value

    @property
    def sigma8(self):
        """
        The present day value of ``sigma_r(r=8 Mpc/h)``, used to normalize
        the power spectrum, which is proportional to the square of this value.

        The power spectrum can re-normalized by setting a different
        value for this parameter
        """
        return self._sigma8

    @sigma8.setter
    def sigma8(self, value):
        """
        Set the sigma8 value and normalize the power spectrum to the new value
        """
        self._norm *= (value / self._sigma8)**2
        self._sigma8 = value

    def __call__(self, k):
        """
        Return the linear power spectrum in units of
        :math:`h^{-3} \mathrm{Mpc}^3` at the redshift specified by :attr:`z`.

        The transfer function used to evaluate the power spectrum is
        specified by the ``transfer`` attribute.

        Parameters
        ---------
        k : float, array_like
            the wavenumber in units of :math:`h Mpc^{-1}`

        Returns
        -------
        Pk : float, array_like
            the linear power spectrum evaluated at ``k`` in units of
            :math:`h^{-3} \mathrm{Mpc}^3`
        """
        if self.transfer != "CLASS":
            Pk = k**self.cosmo.n_s * self._transfer(k)**2
        else:
            k = numpy.asarray(k)
            kmax = self.cosmo.P_k_max
            inrange = k < 0.99999*kmax # prevents rounding errors

            # the return array (could be scalar array)
            Pk = numpy.zeros_like(k)

            # k values in and out of valid range
            k_in = k[inrange]; k_out = k[~inrange]

            # use CLASS in range
            Pk[inrange] = k_in**self.cosmo.n_s * self._transfer(k_in)**2

            # use Eisentein-Hu out of range
            if len(k_out):
                analytic_Tk = self._fallback(k_out)
                analytic_Tk *= self._transfer(kmax)/ self._fallback(kmax)
                Pk[~inrange] = k_out**self.cosmo.n_s * analytic_Tk**2

        return self._norm * Pk

    def velocity_dispersion(self, kmin=1e-5, kmax=10.):
        r"""
        The velocity dispersion in units of of :math:`\mathrm{Mpc/h}`.

        This returns :math:`\sigma_v`, defined as

        .. math::

            \sigma_v^2 = \frac{1}{3} \int_a^b \frac{d^3 q}{(2\pi)^3} \frac{P(q)}{q^2}.

        Parameters
        ----------
        kmin : float, optional
            the lower bound for the integral, in units of :math:`\mathrm{Mpc/h}`
        kmax : float, optional
            the upper bound for the integral, in units of :math:`\mathrm{Mpc/h}`
        """
        from scipy.integrate import quad

        def integrand(logq):
            q = numpy.exp(logq)
            return q*self(q)
        sigmasq = quad(integrand, numpy.log(kmin), numpy.log(kmax))[0] / (6*numpy.pi**2)
        return sigmasq**0.5

    def sigma_r(self, r, kmin=1e-5, kmax=1e1):
        r"""
        The mass fluctuation within a sphere of radius ``r``, in
        units of :math:`h^{-1} Mpc`. This returns :math:`\sigma`, where

        .. math::

            \sigma^2 = \int_0^\infty \frac{k^3 P(k)}{2\pi^2} W^2_T(kr) \frac{dk}{k},

        where :math:`W_T(x) = 3/x^3 (\mathrm{sin}x - x\mathrm{cos}x)` is
        a top-hat filter in Fourier space.

        The value of this function with ``r=8`` returns
        :attr:`sigma8`, within numerical precision.

        Parameters
        ----------
        r : float, array_like
            the scale to compute the mass fluctation over, in units of
            :math:`h^{-1} Mpc`
        kmin : float, optional
            the lower bound for the integral, in units of :math:`\mathrm{Mpc/h}`
        kmax : float, optional
            the upper bound for the integral, in units of :math:`\mathrm{Mpc/h}`
        """
        import mcfit
        from scipy.interpolate import InterpolatedUnivariateSpline as spline

        k = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), 1024)
        Pk = self(k)
        R, sigmasq = mcfit.TophatVar(k)(Pk)

        return spline(R, sigmasq)(r)**0.5
