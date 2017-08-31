from classylss.binding import *
import astropy.cosmology as ac
import astropy.units as units
import numpy
from numpy.testing import assert_allclose
import pytest

# define the astropy cosmologies
c_flat = ac.FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.04, Tcmb0=2.7255)
c_open = ac.LambdaCDM(H0=70, Om0=0.3, Ob0=0.04, Ode0=0.65, Tcmb0=2.7255)
c_closed = ac.LambdaCDM(H0=70, Om0=0.3, Ob0=0.04, Ode0=0.85, Tcmb0=2.7255)

@pytest.mark.parametrize("c2, a_max", [
             (c_flat, 1.0), (c_open, 1.0), (c_closed, 1.0),
             (c_flat, 2.0), (c_open, 2.0), (c_closed, 2.0),
            ])
def test_against_astropy(c2, a_max):
    if a_max > 1.0:
        z = numpy.array([1.0, 0., -0.25])
    else:
        z = numpy.array([1.0, 0.])

    # classylss
    pars = {'h':c2.h, 'Omega_cdm':c2.Odm0, 'Omega_b':c2.Ob0, 'Omega_k':c2.Ok0, 'T_cmb':c2.Tcmb0.value, 'a_max':a_max}
    cosmo = ClassEngine(pars)
    ba = Background(cosmo)

    # age
    this = ba.time(z)
    that = c2.age(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # comoving distance
    this = ba.comoving_distance(z) # in Mpc/h
    that = c2.comoving_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # angular diameter distance
    this = ba.angular_diameter_distance(z) # in Mpc/h
    that = c2.angular_diameter_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)


    # comoving transverse distance
    this = ba.comoving_transverse_distance(z) # in Mpc/h
    that = c2.comoving_transverse_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # critical density
    this = ba.rho_crit(z) # in 1e10 M_sun h^2 / Mpc^3
    that = c2.critical_density(z)
    that = that.to((units.Msun)/(units.Mpc)**3).value / c2.h**2 / 1e10
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # efunc
    this = ba.efunc(z)
    that = c2.efunc(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # luminosity distance
    this = ba.luminosity_distance(z) # in Mpc/h
    that = c2.luminosity_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # baryon omega
    this1 = ba.Omega_b(z)
    this2 = ba.rho_b(z) / ba.rho_crit(z)
    that = c2.Ob(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # lambda omega
    this1 = ba.rho_lambda(z) / ba.rho_crit(z)
    this2 = ba.Omega_lambda(z)
    that = c2.Ode(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # cdm omega
    this1 = ba.rho_cdm(z) / ba.rho_crit(z)
    this2 = ba.Omega_cdm(z)
    that = c2.Odm(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # photons omega
    this1 = ba.rho_g(z) / ba.rho_crit(z)
    this2 = ba.Omega_g(z)
    that = c2.Ogamma(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # total radiation omega
    this1 = (ba.rho_g(z) + ba.rho_ur(z)) / ba.rho_crit(z)
    this2 = ba.Omega_g(z) + ba.Omega_ur(z)
    that = c2.Ogamma(z) + c2.Onu(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # curvature omega
    this1 = ba.rho_k(z) / ba.rho_crit(z)
    this2 = ba.Omega_k(z)
    that = c2.Ok(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # matter omega
    this1 = ba.rho_m(z) / ba.rho_crit(z)
    this2 = ba.Omega_m(z)
    that = c2.Om(z)
    assert_allclose(this1, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
    assert_allclose(this2, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # T cmb
    this = ba.T_cmb(z)
    that = c2.Tcmb(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
