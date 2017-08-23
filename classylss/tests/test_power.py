from classylss.power import LinearPower
from classylss.cosmology import Cosmology
import numpy
from numpy.testing import assert_allclose

def test_linear():

    # initialize the power
    c = Cosmology()
    c.sigma8 = 0.82
    P = LinearPower(c, z=0, transfer='CLASS')

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

    Pk1 = P(k)
    Pk2 = P2(k)
    Pk3 = P3(k)

    # desired redshift
    z = 0

    # linear power spectrum in [Mpc/h]^3
    Plin = power.linear(k, z, verbose=True, cosmo=Planck15)

    # nonlinear power spectrum in [Mpc/h]^3
    Pnl = power.nonlinear(k, z, verbose=True, cosmo=Planck15)

    # Zeldovich power spectrum in [Mpc/h]^3
    Pzel = power.zeldovich(k, z, verbose=True, cosmo=Planck15)
