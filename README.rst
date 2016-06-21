classy_lss
----------
A python binding of the CMB Boltzmann code `CLASS`_, designed for large-scale structure calculations.

The package mainly uses CLASS in order to compute linear power spectra, but the python binding also supports computing nonlinear power 
spectra and the Cls spectra. 

Some of the package features:

- underyling calculations done in C++, with SWIG used to produce the Python interface
- a parameter interface with the `astropy.cosmology`_ module
- linear correlation function calculations using `FFTLog`_
- power spectra and correlation functions in the Zel'dovich approximation

.. _`CLASS` : http://class-code.net
.. _`astropy.cosmology` : http://docs.astropy.org/en/latest/cosmology/index.html
.. _`FFTLog` : http://casa.colorado.edu/~ajsh/FFTLog/

Dependencies
------------

The required external dependencies are: 

- FFTW
- astropy
- numpy
- requests

The latest version of CLASS will automatically be downloaded and compiled, and is thus, not an external dependency for the user.

Installation
------------

The package can be cloned from github as

.. code:: bash

    git clone https://github.com/nickhand/classy_lss.git
    cd classy_lss
    
If the ``fftw3`` library is installed in a non-standard location, the user should specify the
location using environment variables:

.. code:: bash

    export FFTW_INC=/path/to/fftw/include
    export FFTW_DIR=/path/to/fftw/lib

and then the ``classy_lss`` package can be installed using

.. code:: bash
    
    export CFLAGS='-std=c++11 -fopenmp'; pip install .
    
This procedure has been tested on Mac and Linux machines.

To verify that the installation has succeeded, run:

.. code-block:: python

    from classy_lss import gcl
    
Examples
--------

To compute power spectra for the Planck 2015 cosmology:

.. code-block:: python

    from astropy.cosmology import Planck15
    from classy_lss import power
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
    
and similarly, correlation functions can be computed: 

.. code-block:: python
    
    from classy_lss import correlation
    
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
    
All of the above functions accept a ``class_kwargs`` keyword, which allows the user
to pass any valid CLASS parameter to the CLASS code. The ``class_kwargs`` parameter is a dictionary 
that will be passed to the ``ClassEngine`` instance, which is responsible for running CLASS. 

    
