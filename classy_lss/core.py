from . import ClassParams, LinearPS, Cosmology, ClassEngine
from . import transfers

import numpy

__all__ = ['linear_power', 'nonlinear_power']

def linear_power(k, z, transfer=transfers.CLASS, cosmo=None, verbose=False, class_kws={}):
    """
    Compute the linear power spectrum P_lin(k, z) using the specified 
    transfer function
    
    Parameters
    ----------
    k : array_like
        the wavenumbers in units of `h/Mpc`
    z : float
        the redshift to compute the power at
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
    Plin : array_like
        the linear power spectrum in units of (Mpc/h)^3
    """
    from astropy.cosmology import default_cosmology, FLRW
    if cosmo is None:
        cosmo = default_cosmology.get()
    
    if not isinstance(cosmo, FLRW):
        raise TypeError("cosmology should be a astropy.cosmology.FLRW subclass")
    if isinstance(k, list):
        k = numpy.array(k)
    
    # default CLASS parameters
    class_kws.setdefault('P_k_max_h/Mpc',  20.)
    class_kws.setdefault('z_max_pk', 2.0)
    
    # astropy to CLASS params
    pars = ClassParams.from_astropy(cosmo, extra=class_kws)
    
    # do the work
    class_cosmo = Cosmology(pars, transfer, verbose)
    Plin = LinearPS(class_cosmo, z)
    
    return Plin(k)    
    
def nonlinear_power(k, z, cosmo=None, verbose=False, class_kws={}):
    """
    Compute the nonlinear power spectrum P_nl(k, z) using HALOFIT
    
    Parameters
    ----------
    k : array_like
        the wavenumbers in units of `h/Mpc`
    z : float
        the redshift to compute the power at
    cosmo : astropy.cosmology
        the astropy cosmology; defaults to ``default_cosmology()``
    verbose : bool, optional
        whether to be verbose when running CLASS; default is `False`
    class_kws : dict, optional
        extra parameter keywords to use when running CLASS
    
    Returns
    -------
    Plin : array_like
        the linear power spectrum in units of (Mpc/h)^3
    """
    from astropy.cosmology import default_cosmology, FLRW
    if cosmo is None:
        cosmo = default_cosmology.get()
    
    if not isinstance(cosmo, FLRW):
        raise TypeError("cosmology should be a astropy.cosmology.FLRW subclass")
    if isinstance(k, list):
        k = numpy.array(k)
    
    # set the defaul class keywords
    class_kws.setdefault('non linear', 'halofit')
    class_kws.setdefault('P_k_max_h/Mpc',  20.)
    class_kws.setdefault('z_max_pk', 2.0)
    
    # astropy to CLASS params
    pars = ClassParams.from_astropy(cosmo, extra=class_kws)
    
    # make sure we are computing mPk
    if 'output' not in pars['output']:
        pars['output'] = "mPk"
    elif 'mPk' not in pars['output']:
        pars['output'] += ", mPk"
    
    # do the work
    class_cosmo = ClassEngine(verbose)
    with class_cosmo.compute(pars):
        Pnl = class_cosmo.GetPknl(k, z)
    
    return Pnl   
    
    
    