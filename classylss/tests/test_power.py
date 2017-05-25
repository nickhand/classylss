
def test_power():

    from astropy.cosmology import Planck15
    from classylss import power
    import numpy
    
    # desired wavenumbers (in h/Mpc)
    k = numpy.logspace(-3, 0, 500)

    # desired redshift
    z = 0

    # linear power spectrum in [Mpc/h]^3
    Plin = power.linear(k, z, verbose=True, cosmo=Planck15)

    # nonlinear power spectrum in [Mpc/h]^3
    Pnl = power.nonlinear(k, z, verbose=True, cosmo=Planck15)

    # Zeldovich power spectrum in [Mpc/h]^3
    Pzel = power.zeldovich(k, z, verbose=True, cosmo=Planck15)