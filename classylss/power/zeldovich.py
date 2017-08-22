import mcfit
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.integrate import quad

from .linear import LinearPower

NUM_PTS = 1024
KMIN = 1e-5
KMAX = 1e1
NMAX = 15

def isiterable(obj):
    """Returns `True` if the given object is iterable."""

    try:
        iter(obj)
        return True
    except TypeError:
        return False

def vectorize_if_needed(func, *x):
    """ Helper function to vectorize functions on array inputs"""
    if any(map(isiterable, x)):
        return numpy.vectorize(func)(*x)
    else:
        return func(*x)

class ZeldovichPower(object):
    """
    The matter power spectrum in the Zel'dovich approximation.
    """

    def __init__(self, cosmo, z, transfer='CLASS'):

        # initialize the linear power
        self.Plin = LinearPower(cosmo, z, transfer='CLASS')

        self.cosmo = cosmo
        self._sigma8 = self.cosmo.sigma8
        self.z = z

        # use low-k approx below this value
        self._k0_low = 5e-3

    def _setup(self):

        # set up the k-grid for integrals
        k = numpy.logspace(numpy.log10(KMIN), numpy.log10(KMAX), NUM_PTS)
        Pk = self.Plin(k)

        # compute the I0, I1 integrals we need
        self._r, I0 = ZeldovichJ0(k)(Pk, extrap=False)
        _, I1 = ZeldovichJ1(k)(Pk, extrap=False)

        # compute the X(r), Y(r) integrals we need
        self._sigmasq = self.Plin.velocity_dispersion(kmin=KMIN, kmax=KMAX)**2
        self._X = -2.*I1 + 2 * self._sigmasq
        self._Y = -2.*I0 + 6.*I1

        # needed for the low-k approx
        self._Q3 = quad(lambda q: (self.Plin(q)/q)**2, 1e-6, 100.)[0]

    @property
    def z(self):
        """
        The redshift of the power spectrum
        """
        return self._z

    @z.setter
    def z(self, value):
        self._z = value
        self.Plin.z = value
        self._setup()

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
        self._sigma8 = value
        self.Plin.sigma8 = value
        self._setup()

    def _low_k_approx(self, k):
        r"""
        Return the low-k approximation of the Zel'dovich power. This computes:

        .. math::
            P(k) = (1 - k^2 \sigma_v^2 + 1/2 k^4 \sigma_v^4) P_L(k) + 0.5 Q_3(k),

        where :math:`Q_3(k)` is defined as

        .. math::
            Q_3(k) = \frac{k^4}{10\pi^2} \int dq \frac{P^2_L(q)}{q^2}.
        """
        Q3 = 1./(10.*numpy.pi**2) * k**4 * self._Q3

        Plin = self.Plin(k)
        term1 = (1 - k**2 * self._sigmasq + 0.5 * k**4 * self._sigmasq**2) * Plin
        term2 = 0.5 * Q3
        return term1 + term2

    def __call__(self, k):
        """
        Return the Zel'dovich power in :math:`h^{-3} \mathrm{Mpc}^3 at
        :attr:`z` and ``k``, where ``k`` is in units of
        :math:`h \mathrm{Mpc}^{-1}`.

        Parameters
        ----------
        k : float, array_like
            the wavenumbers to evaluate the power at
        """
        def Pzel_at_k(ki):

            # return the low-k approximation
            if ki < self._k0_low:
                return self._low_k_approx(ki)

            # do the full integral
            Pzel = numpy.zeros(len(self._r))
            for n in range(0, NMAX+1):

                I = ZeldovichPowerIntegral(self._r, n)
                if n > 0:
                    f = self._Y**n * numpy.exp(-0.5*ki**2 * (self._X + self._Y))
                else:
                    f = numpy.exp(-0.5*ki**2 * (self._X + self._Y)) - numpy.exp(-ki**2*self._sigmasq)

                kk, this_Pzel = I(f, extrap=False)
                Pzel += kk**n * this_Pzel

            return spline(kk, Pzel)(ki)

        return vectorize_if_needed(Pzel_at_k, k)

class ZeldovichJ0(mcfit.mcfit):
    """
    An integral over :math:`j_0` needed to compute the Zeldovich power. The
    integral is given by:

    .. math::

        I_0(r) = \int \frac{dk}{2\pi^2} P_L(k) j_0(kr).
    """
    def __init__(self, k, N=None):
        self.l = 0
        UK = mcfit.kernels.Mellin_SphericalBesselJ(self.l)
        prefac = k
        postfac = 1 / (2*numpy.pi)**1.5
        mcfit.mcfit.__init__(self, k, UK, q=1.0, N=N, lowring=False, prefac=prefac, postfac=postfac)

class ZeldovichJ1(mcfit.mcfit):
    """
    An integral over :math:`j_1` needed to compute the Zeldovich power. The
    integral is given by:

    .. math::

        I_1(r) = \int \frac{dk}{2\pi^2} P_L(k) \frac{j_1(kr)}{kr}.
    """
    def __init__(self, k, N=None):
        self.l = 1
        UK = mcfit.kernels.Mellin_SphericalBesselJ(self.l)
        prefac = 1.0
        postfac = 1 / (2*numpy.pi)**1.5
        mcfit.mcfit.__init__(self, k, UK, q=0., N=N, lowring=False, prefac=prefac, postfac=postfac)
        self.postfac /= self.y


class ZeldovichPowerIntegral(mcfit.mcfit):
    """
    The integral needed to evaluate the density auto spectrum in the
    Zel'dovich approximation.

    This evaluates:

    .. math::
        I(k, n) =
    """
    def __init__(self, r, n, N=None):
        self.n = n
        UK = mcfit.kernels.Mellin_SphericalBesselJ(self.n)

        prefac = r**(3-n)
        postfac = (2*numpy.pi)**1.5
        mcfit.mcfit.__init__(self, r, UK, q=1.5, N=N, lowring=False, prefac=prefac, postfac=postfac)
