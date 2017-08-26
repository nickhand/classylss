from classylss.binding import *
import astropy.cosmology as ac
import astropy.units as units
import numpy
from numpy.testing import assert_allclose

def test_flat():

    # astropy
    c2 = ac.FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.04, Tcmb0=2.7255)

    # classylss
    pars = {'h':0.7, 'Omega_cdm':0.26, 'Omega_b':0.04, 'Omega_k':0., 'T_cmb':2.7255, 'a_max':2.0}
    cosmo = ClassEngine(pars)
    ba = Background(cosmo)

    # age
    z = numpy.array([1.0, 0., -0.25])
    this = ba.time(z)
    that = c2.age(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # angular diameter distance
    z = numpy.array([1.0, -0.25])
    this = ba.angular_diameter_distance(z) # in Mpc/h
    that = c2.angular_diameter_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # comoving distance
    z = numpy.array([1.0, -0.25])
    this = ba.comoving_distance(z) # in Mpc/h
    that = c2.comoving_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # comoving transverse distance
    z = numpy.array([1.0, -0.25])
    this = ba.comoving_transverse_distance(z) # in Mpc/h
    that = c2.comoving_transverse_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # critical density
    z = numpy.array([1.0, 0., -0.25])
    this = ba.rho_crit(z) # in 1e10 M_sun h^2 / Mpc^3
    that = c2.critical_density(z)
    that = that.to((units.Msun)/(units.Mpc)**3).value / c2.h**2 / 1e10
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # efunc
    z = numpy.array([1.0, 0., -0.25])
    this = ba.efunc(z)
    that = c2.efunc(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # luminosity distance
    z = numpy.array([1.0, -0.25])
    this = ba.luminosity_distance(z) # in Mpc/h
    that = c2.luminosity_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # baryon omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_b(z)
    that = c2.Ob(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # lambda omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_lambda(z)
    that = c2.Ode(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # cdm omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_cdm(z)
    that = c2.Odm(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # photons omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_g(z)
    that = c2.Ogamma(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # total radiation omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_g(z) + ba.Omega_ur(z)
    that = c2.Ogamma(z) + c2.Onu(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # curvature omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_k(z)
    that = c2.Ok(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # matter omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_m(z)
    that = c2.Om(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # T cmb
    z = numpy.array([1.0, 0., -0.25])
    this = ba.T_cmb(z)
    that = c2.Tcmb(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)


def test_nonflat():

    # astropy
    c2 = ac.LambdaCDM(H0=70, Om0=0.3, Ob0=0.04, Ode0=0.80, Tcmb0=2.7255)

    # classylss
    pars = {'h':0.7, 'Omega_cdm':0.26, 'Omega_b':0.04, 'Omega_k':c2.Ok0, 'T_cmb':2.7255, 'a_max':2.0}
    cosmo = ClassEngine(pars)
    ba = Background(cosmo)

    # age
    z = numpy.array([1.0, 0., -0.25])
    this = ba.time(z)
    that = c2.age(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # angular diameter distance
    z = numpy.array([1.0, -0.25])
    this = ba.angular_diameter_distance(z) # in Mpc/h
    that = c2.angular_diameter_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # comoving distance
    z = numpy.array([1.0, -0.25])
    this = ba.comoving_distance(z) # in Mpc/h
    that = c2.comoving_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # comoving transverse distance
    z = numpy.array([1.0, -0.25])
    this = ba.comoving_transverse_distance(z) # in Mpc/h
    that = c2.comoving_transverse_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # critical density
    z = numpy.array([1.0, 0., -0.25])
    this = ba.rho_crit(z) # in 1e10 M_sun h^2 / Mpc^3
    that = c2.critical_density(z)
    that = that.to((units.Msun)/(units.Mpc)**3).value / c2.h**2 / 1e10
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # efunc
    z = numpy.array([1.0, 0., -0.25])
    this = ba.efunc(z)
    that = c2.efunc(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # luminosity distance
    z = numpy.array([1.0, -0.25])
    this = ba.luminosity_distance(z) # in Mpc/h
    that = c2.luminosity_distance(z).value * c2.h
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # baryon omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_b(z)
    that = c2.Ob(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # lambda omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_lambda(z)
    that = c2.Ode(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # cdm omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_cdm(z)
    that = c2.Odm(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # photons omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_g(z)
    that = c2.Ogamma(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # total radiation omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_g(z) + ba.Omega_ur(z)
    that = c2.Ogamma(z) + c2.Onu(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # curvature omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_k(z)
    that = c2.Ok(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # matter omega
    z = numpy.array([1.0, 0., -0.25])
    this = ba.Omega_m(z)
    that = c2.Om(z)
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)

    # T cmb
    z = numpy.array([1.0, 0., -0.25])
    this = ba.T_cmb(z)
    that = c2.Tcmb(z).value
    assert_allclose(this, that, rtol=1e-3, atol=1e-3, err_msg='z = %s' %z)
