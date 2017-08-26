from classylss.power import LinearPower, HalofitPower, ZeldovichPower
from classylss.cosmology import Cosmology
import numpy
from numpy.testing import assert_allclose
import pytest

def test_large_scales():

    c = Cosmology()
    k = numpy.logspace(-5, -2, 100)

    # linear power
    Plin = LinearPower(c, z=0)

    # nonlinear power
    Pnl = HalofitPower(c, z=0)

    # zeldovich power
    Pzel = ZeldovichPower(c, z=0)

    assert_allclose(Plin(k), Pnl(k), rtol=1e-2)
    assert_allclose(Plin(k), Pzel(k), rtol=1e-2)

def test_linear():

    # initialize the power
    c = Cosmology()
    c.sigma8 = 0.82
    P = LinearPower(c, z=0, transfer='CLASS')

    # check velocity dispersion
    assert_allclose(P.velocity_dispersion(), 5.898, rtol=1e-3)

    # test sigma8
    assert_allclose(P.sigma_r(8.), c.sigma8, rtol=1e-5)

    # change sigma8
    P.sigma8 = c.sigma8 = 0.80
    assert_allclose(P.sigma_r(8.), P.sigma8, rtol=1e-5)

    # change redshift and test sigma8(z)
    P.z = 0.55
    assert_allclose(P.sigma_r(8.), c.sigma8_z(P.z), rtol=1e-5)

    # desired wavenumbers (in h/Mpc)
    k = numpy.logspace(-3, 2, 500)

    # initialize EH power
    P1 = LinearPower(c, z=0., transfer="CLASS")
    P2 = LinearPower(c, z=0., transfer='EisensteinHu')
    P3 = LinearPower(c, z=0., transfer='NoWiggleEisensteinHu')

    # check different transfers (very roughly)
    Pk1 = P1(k)
    Pk2 = P2(k)
    Pk3 = P3(k)
    assert_allclose(Pk1 / Pk1.max(), Pk2 / Pk2.max(), rtol=0.1)
    assert_allclose(Pk1 / Pk1.max(), Pk3 / Pk3.max(), rtol=0.1)

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
