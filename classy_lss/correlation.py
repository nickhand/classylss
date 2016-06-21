from . import power, tools
import numpy

NUM_PTS = 1024

def linear(r, z, kmin=1e-5, kmax=10., smoothing=0., **power_kws):
    """
    Compute the linear correlation function using the specified 
    transfer function
    
    Parameters
    ----------
    r : array_like
        the separation in units of `Mpc/h`
    z : float
        the redshift to compute the power at
    kmin : float, optional
        minimum wavenumber in `h/Mpc` to integrate over when computing CF
    kmax : float, optional
        maximum wavenumber in `h/Mpc` to integrate over when computing CF
    smoothing : float, optional
        smoothing scale in `Mpc/h`
    **power_kws : 
        keywords passed to :func:`power.linear`
            transfer : int, optional
                available transfer functions defined in `classy_lss.transfers`;
                the default value runs CLASS
            cosmo : astropy.cosmology, optional
                the astropy cosmology; defaults to ``default_cosmology()``
            verbose : bool, optional
                whether to be verbose when running CLASS; default is `False`
            class_kws : dict, optional
                extra parameter keywords to use when running CLASS
    
    Returns
    -------
    cf : array_like
        the linear correlation function (unitless)
    """
    k = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), NUM_PTS)
    Pk = power.linear(k, z, **power_kws)
    
    return tools.pk_to_xi(0, k, Pk, r, smoothing);
    
def nonlinear(r, z, kmin=1e-5, kmax=10., smoothing=0., **power_kws):
    """
    Compute the nonlinear correlation function using HALOFIT

    Parameters
    ----------
    r : array_like
        the separation in units of `Mpc/h`
    z : float
        the redshift to compute the power at
    kmin : float, optional
        minimum wavenumber in `h/Mpc` to integrate over when computing CF
    kmax : float, optional
        maximum wavenumber in `h/Mpc` to integrate over when computing CF
    smoothing : float, optional
        smoothing scale in `Mpc/h`
    **power_kws : 
        keywords passed to :func:`power.nonlinear`
            cosmo : astropy.cosmology
                the astropy cosmology; defaults to ``default_cosmology()``
            verbose : bool, optional
                whether to be verbose when running CLASS; default is `False`
            class_kws : dict, optional
                extra parameter keywords to use when running CLASS

    Returns
    -------
    cf : array_like
        the nonlinear correlation function (unitless)
    """
    k = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), NUM_PTS)
    Pk = power.nonlinear(k, z, **power_kws)
    
    return tools.pk_to_xi(0, k, Pk, r, smoothing);

def zeldovich(r, z, kmin=1e-5, kmax=10., smoothing=0., **power_kws):
    """
    Compute the Zel'dovich correlation function using the specified
    transfer function for the Zel'dovich power spectrum

    Parameters
    ----------
    r : array_like
        the separation in units of `Mpc/h`
    z : float
        the redshift to compute the power at
    kmin : float, optional
        minimum wavenumber in `h/Mpc` to integrate over when computing CF
    kmax : float, optional
        maximum wavenumber in `h/Mpc` to integrate over when computing CF
    smoothing : float, optional
        smoothing scale in `Mpc/h`
    **power_kws : 
        keywords passed to :func:`power.zeldovich`
            transfer : int, optional
                available transfer functions defined in `classy_lss.transfers`;
                the default value runs CLASS
            cosmo : astropy.cosmology, optional
                the astropy cosmology; defaults to ``default_cosmology()``
            verbose : bool, optional
                whether to be verbose when running CLASS; default is `False`
            class_kws : dict, optional
                extra parameter keywords to use when running CLASS

    Returns
    -------
    cf : array_like
        the Zel'dovich correlation function (unitless)
    """
    k = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), NUM_PTS)
    Pk = power.zeldovich(k, z, **power_kws)
    
    return tools.pk_to_xi(0, k, Pk, r, smoothing);
