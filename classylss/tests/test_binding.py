from classylss.binding import *
import pytest

@pytest.mark.parametrize("a_max", [1.0, 2.0])
def test_ba(a_max):
    cosmo = ClassEngine({'a_max':a_max, 'N_ncdm': 1, 'm_ncdm':[0.06]})
    ba = Background(cosmo)
    assert ba.hubble_function(0.1).shape == ()
    assert ba.hubble_function([0.1]).shape == (1,)
    assert ba.hubble_function([0.2, 0.3]).shape == (2,)
    assert ba.hubble_function([[0.2, 0.3, 0.4]]).shape == (1, 3,)
    ba.time([[0.2, 0.3, 0.4]])
    ba.comoving_distance([[0.2, 0.3, 0.4]])
    ba.angular_diameter_distance([[0.2, 0.3, 0.4]])
    ba.luminosity_distance([[0.2, 0.3, 0.4]])
    ba.hubble_function([[0.2, 0.3, 0.4]])
    ba.hubble_function_prime([[0.2, 0.3, 0.4]])
    ba.scale_independent_growth_factor([[0.2, 0.3, 0.4]])
    ba.scale_independent_growth_rate([[0.2, 0.3, 0.4]])

    assert ba.Omega_pncdm(0, 0) == ba.Omega0_pncdm[0]

@pytest.mark.parametrize("a_max", [1.0, 2.0])
def test_sp(a_max):
    cosmo = ClassEngine({'a_max':a_max, 'output': 'dTk vTk mPk', 'P_k_max_h/Mpc' : 20., "z_max_pk" : 100.0})
    sp = Spectra(cosmo)
    pk = sp.get_pk(0.1, 0.1)
    pk = sp.get_pk([0.1], 0.1)
    pk = sp.get_pk([0.1], [0.1])

    t = sp.get_transfer(0.0)
    t = sp.get_transfer(2.0)
    t = sp.get_transfer(1./(0.99*a_max)-1)

    s8_z = sp.sigma8_z([0., 1.0])

@pytest.mark.parametrize("a_max", [1.0, 2.0])
def test_thermo(a_max):
    cosmo = ClassEngine({'a_max':a_max, 'output': 'dTk vTk mPk', 'P_k_max_h/Mpc' : 20., "z_max_pk" : 100.0})
    th = Thermo(cosmo)

    z_d = th.z_drag
    rs_d = th.rs_drag
    tau_reio = th.tau_reio
    z_reio = th.z_reio
    z_rec = th.z_rec
    rs_res = th.rs_rec
    theta_s = th.theta_s

@pytest.mark.parametrize("a_max", [1.0, 2.0])
def test_primordial(a_max):

    cosmo = ClassEngine({'a_max':a_max, 'output': 'dTk vTk mPk', 'P_k_max_h/Mpc' : 20., "z_max_pk" : 100.0})
    pm = Primordial(cosmo)
    Pk = pm.get_pkprim([0., 0.1, 0.2])

    pr = pm.get_primordial()
