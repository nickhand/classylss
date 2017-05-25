def test_corr():
    
    from astropy.cosmology import Planck15
    from classylss import correlation
    import numpy
    
    # desired separation (in Mpc/h)
    r = numpy.logspace(0, numpy.log10(150), 500)

    # desired redshift
    z = 0

    # linear 2PCF
    cf_lin = correlation.linear(r, z, verbose=True, cosmo=Planck15)

    # nonlinear 2PCF
    cf_nl = correlation.nonlinear(r, z, verbose=True, cosmo=Planck15)

    # Zeldovich power spectrum in [Mpc/h]^3
    cf_zel = correlation.zeldovich(r, z, smoothing=1.0, verbose=True, cosmo=Planck15)