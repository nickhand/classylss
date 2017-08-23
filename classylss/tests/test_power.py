from classylss.power import LinearPower, HalofitPower, ZeldovichPower
from classylss.cosmology import Cosmology
import numpy
from numpy.testing import assert_allclose
import pytest

def test_linear():

    # initialize the power
    c = Cosmology()
    c.sigma8 = 0.82
    P = LinearPower(c, z=0, transfer='CLASS')

    # check velocity dispersion
    assert_allclose(P.velocity_dispersion(), 6.034, rtol=1e-3)

    # test sigma8
    assert_allclose(P.sigma_r(8.), c.sigma8, rtol=1e-3)

    # change sigma8
    c.sigma8 = 0.80
    assert_allclose(P.sigma_r(8.), c.sigma8, rtol=1e-3)

    # change redshift and test sigma8(z)
    P.z = 0.55
    assert_allclose(P.sigma_r(8.), c.sigma8_z(P.z), rtol=1e-3)

    # desired wavenumbers (in h/Mpc)
    k = numpy.logspace(-3, 2, 500)

    # initialize EH power
    P2 = LinearPower(c, z=P.z, transfer='EisensteinHu')
    P3 = LinearPower(c, z=P.z, transfer='NoWiggleEisensteinHu')

    # check different transfers
    Pk1 = P(k)
    Pk2 = P2(k)
    Pk3 = P3(k)
    assert_allclose(P1k, P2k)
    assert_allclose(P1k, P3k)

    # also try scalar
    Pk = P(0.1)

def test_halofit():

    # initialize the power
    c = Cosmology()
    c.sigma8 = 0.82
    P = HalofitPower(c, z=0)

    # k is out of range
    with pytest.raises(ValueError):
        Pk = P(2*c.P_k_max)

    # compute for scalar
    Pk = P(0.1)

    # compute for array
    k = numpy.logspace(-3, numpy.log10(c.P_k_max), 500)
    Pk = P(k)

def test_zeldovich():

    # initialize the power
    c = Cosmology()
    c.sigma8 = 0.82
    P = ZeldovichPower(c, z=0)

    # compute for scalar
    Pk = P(0.1)

    # compute for array
    k = numpy.logspace(-3, numpy.log10(c.P_k_max), 500)
    Pk = P(k)
