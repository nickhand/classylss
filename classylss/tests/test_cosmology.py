from classylss.cosmology import Cosmology
from numpy.testing import assert_allclose, assert_array_equal
from numpy import meshgrid

def test_load_precision():
    from classylss import load_precision

    p = load_precision('pk_ref.pre')

    c = Cosmology(gauge='synchronous', tol_background_integration=1e-5, **p)

    assert_allclose(c.Omega_cdm(0), c.Omega0_cdm)

def test_cosmology_sane():
    c = Cosmology(gauge='synchronous')

    assert_allclose(c.Omega_cdm(0), c.Omega0_cdm)
    assert_allclose(c.Omega_g(0), c.Omega0_g)
    assert_allclose(c.Omega_b(0), c.Omega0_b)
    assert_allclose(c.Omega_ncdm(0), c.Omega0_ncdm)
    assert_allclose(c.Omega_ur(0), c.Omega0_ur)
    assert_allclose(c.Omega_ncdm(0), c.Omega0_ncdm_tot)

    assert_allclose(c.Omega_pncdm(0), c.Omega0_pncdm)
    assert_allclose(c.Omega_m(0), c.Omega0_m)
    assert_allclose(c.Omega_r(0), c.Omega0_r)

    # total density in 10e10 Msun/h unit
    assert_allclose(c.rho_tot(0), 27.754999)

    # comoving distance to z=1.0 in Mpc/h unit.
    assert_allclose(c.conformal_distance(1.0), 3396.157391 * c.h)

    # conformal time in Mpc unit.
    assert_allclose(c.tau(1.0), 3396.157391)

    assert_allclose(c.efunc(0), 1.) # hubble in Mpc/h km/s unit
    assert_allclose(c.efunc(0) - c.efunc(1 / 0.9999 - 1),
                    0.0001 * c.efunc_prime(0), rtol=1e-3)

def test_cosmology_density():
    c = Cosmology(gauge='synchronous')
    z = [0, 1, 2, 5, 9, 99]
    assert_allclose(c.rho_cdm(z), c.Omega_cdm(z) * c.rho_tot(z))
    assert_allclose(c.rho_g(z), c.Omega_g(z) * c.rho_tot(z))
    assert_allclose(c.rho_ncdm(z), c.Omega_ncdm(z) * c.rho_tot(z))
    assert_allclose(c.rho_b(z), c.Omega_b(z) * c.rho_tot(z))
    assert_allclose(c.rho_m(z), c.Omega_m(z) * c.rho_tot(z))
    assert_allclose(c.rho_r(z), c.Omega_r(z) * c.rho_tot(z))
    assert_allclose(c.rho_ur(z), c.Omega_ur(z) * c.rho_tot(z))

def test_cosmology_vect():
    c = Cosmology(gauge='synchronous')

    assert_allclose(c.Omega_cdm([0]), c.Omega0_cdm)

    assert_array_equal(c.Omega_cdm([]).shape, [0])
    assert_array_equal(c.Omega_cdm([0]).shape, [1])
    assert_array_equal(c.Omega_cdm([[0]]).shape, [1, 1])

    assert_array_equal(c.rho_k([[0]]).shape, [1, 1])

    k, z = meshgrid([0, 1], [0.01, 0.05, 0.1, 0.5], sparse=True, indexing='ij')

    pk = c.get_pk(z=z, k=k)
    assert_array_equal(pk.shape, [2, 4])

def test_cosmology_a_max():
    c = Cosmology(gauge='synchronous', a_max=2.0)
    print(c.parameter_file)
    assert c.a_max == 2.0
    t = c.Omega_m(-0.1)
    t = c.efunc(-0.1)
    t = c.scale_independent_growth_factor(-0.1)

    #FIXME: transfer doesn't work because thermal dynamics doesn't go beyond z=0.
    #t = c.get_transfer(z=-0.1)

def test_cosmology_transfer():
    c = Cosmology()
    t = c.get_transfer(z=0)
    assert 'h_prime' in t.dtype.names
    assert 'k' in t.dtype.names
    assert 'd_cdm' in t.dtype.names

def test_cosmology_get_pk():
    c = Cosmology()
    p = c.get_pk(z=0, k=0.1)
    p1 = c.Spectra.get_pk(z=0, k=0.1)

    # ensure the dro did use get_pk of Spectra rather than that from Primordial
    assert_allclose(p, p1)

def test_cosmology_astropy():
    from astropy.cosmology import Planck15
    c = Cosmology.from_astropy(Planck15)

def test_cosmology_dir():
    c = Cosmology()
    d = dir(c)
    assert "Background" in d
    assert "Spectra" in d
    assert "Omega0_m" in d

def test_cosmology_pickle():
    import pickle
    c = Cosmology()
    s = pickle.dumps(c)
    c1 = pickle.loads(s)
    assert c1.parameter_file == c.parameter_file

def test_cosmology_clone():
    c = Cosmology(gauge='synchronous')

    c1 = Cosmology(gauge='newtonian')
    assert 'newtonian' in c1.parameter_file

    c2 = Cosmology(P_k_max=1.01234567)
    assert '1.01234567' in c2.parameter_file

def test_cosmology_astropy():
    c = Cosmology(gauge='synchronous', m_ncdm=[0.06])

    assert_allclose(c.Odm(0), c.Odm0)
    assert_allclose(c.Ogamma(0), c.Ogamma0)
    assert_allclose(c.Ob(0), c.Ob0)
    assert_allclose(c.Onu(0), c.Onu0)
    assert_allclose(c.Ok(0), c.Ok0)
    assert_allclose(c.Ode(0), c.Ode0)
    assert_array_equal(c.has_massive_nu, True)
