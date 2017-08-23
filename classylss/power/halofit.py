import numpy

class HalofitPower(object):
    """
    Nonlinear power spectrum computed using HaloFit via CLASS.

    Parameters
    ----------
    cosmo : :class:`Cosmology`
        the cosmology instance
    z : float
        the redshift of the power spectrum
    """

    def __init__(self, cosmo, z):

        # make sure nonlinear is enabled in CLASS
        if not cosmo.nonlinear:
            cosmo = cosmo.clone(**{'non linear':'halofit'})
        self.cosmo = cosmo
        self.z = z

    def __call__(self, k):
        r"""
        Return the power in units of :math:`h^{-3} \mathrm{Mpc}^3`.

        Parameters
        ----------
        k : float, array_like
            the wavenumbers in units of :math:`h \mathrm{Mpc}^{-1}`

        Returns
        -------
        Pk : float, array_like
            the linear power spectrum evaluated at ``k`` in units of
            :math:`h^{-3} \mathrm{Mpc}^3`
        """
        k = numpy.asarray(k)
        if k.max() > self.cosmo.P_k_max:
            msg = "results can only be computed up to k=%.2e h/Mpc; " %self.cosmo.P_k_max
            msg += "try increasing the Cosmology parameter 'P_k_max'"
            raise ValueError(msg)

        kmin = self.cosmo.P_k_min
        inrange = k > 1.00001*kmin

        Pk = numpy.zeros_like(k)
        k_in = k[inrange]; k_out = k[~inrange]

        Pk[inrange] = self.cosmo.get_pk(k=k_in, z=self.z)
        Pk[~inrange] = self.cosmo.get_pk(k=k_out, z=self.z)
        return Pk
