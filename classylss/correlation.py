import numpy
import mcfit
from scipy.interpolate import InterpolatedUnivariateSpline
from .power.zeldovich import ZeldovichPower

NUM_PTS = 1024

class CorrelationFunction(object):
    """
    Evaluate the correlation function by Fourier transforming
    a power spectrum object.

    Parameters
    ----------
    power : callable
         a callable power spectrum that returns the power at a given ``k``;
         this should have ``z`` and ``sigma8`` attributes
    """
    def __init__(self, power):

        self.power = power

        if not hasattr(power, 'z'):
            raise AttributeError("input power spectrum must have a ``z`` attribute")
        self._z = power.z

        if not hasattr(power, 'sigma8'):
            raise AttributeError("input power spectrum must have a ``z`` attribute")
        self._sigma8 = power.sigma8

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, value):
        self._z = value
        self.power.z = value

    @property
    def sigma8(self):
        return self._sigma8

    @sigma8.setter
    def sigma8(self, value):
        self._sigma8 = value
        self.power.sigma8 = value

    def __call__(self, r, smoothing=0., kmin=1e-5, kmax=10.):
        """
        Return the correlation function (dimensionless) for separations ``r``

        Parameters
        ----------
        r : float, array_like
            the separation array in units of :math:`h^{-1} \mathrm(Mpc)`
        smoothing  : float, optional
            the std deviation of the Fourier space Gaussian smoothing to apply
            to P(k) before taking the FFT
        kmin : float, optional
            the minimum ``k`` value to compute P(k) for before taking the FFT
        kmax : float, optional
            the maximum ``k`` value to compute P(k) for before taking the FFT
        """
        k = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), NUM_PTS)
        xi = mcfit.P2xi(k)

        # power with smoothing
        Pk = self.power(k)
        Pk *= numpy.exp(-(k*smoothing)**2)

        # only extrap if not zeldovich
        extrap = not isinstance(self.power, ZeldovichPower)
        rr, CF = xi(Pk, extrap=extrap)
        return InterpolatedUnivariateSpline(rr, CF)(r)
