from classylss.cosmology import Cosmology
from numpy.testing import assert_allclose

def test_cosmology_init():
    c = Cosmology(gauge='synchronous')

    assert_allclose(c.Ocdm(0), c.Ocdm0)
    assert_allclose(c.Og(0), c.Og0)
    assert_allclose(c.Ob(0), c.Ob0)
    assert_allclose(c.Oncdm(0), c.Oncdm0)
    assert_allclose(c.Our(0), c.Our0)

    assert_allclose(c.Opncdm(0), c.Opncdm0)
    assert_allclose(c.Om(0), c.Om0) 
    assert_allclose(c.Or(0), c.Or0)

    # total density in 10e10 Msun/h unit
    assert_allclose(c.rho_tot(0), 27.754999)

    # comoving distance to z=1.0 in Mpc/h unit.
    assert_allclose(c.conformal_distance(1.0), 3396.157391 * c.h)

    # conformal time in Mpc unit.
    assert_allclose(c.tau(1.0), 3396.157391)

    assert_allclose(c.efunc(0), 1.) # hubble in Mpc/h km/s unit
    assert_allclose(c.efunc(0) - c.efunc(1 / 0.9999 - 1), 
                    0.0001 * c.efunc_prime(0), rtol=1e-3)

def test_cosmology_a_max():
    c = Cosmology(gauge='synchronous', a_max=2.0)
    print(c.parameter_file)
    assert c.a_max == 2.0
    t = c.Om(-0.1)
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
    assert "Om0" in d

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
