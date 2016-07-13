classylss
----------

.. image:: https://img.shields.io/pypi/v/classylss.svg
   :alt: PyPi
   :target: https://pypi.python.org/pypi/classylss/


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

The required external Python dependencies are: 

- astropy
- numpy

And the necessary compilation tools are: 

- gfortran
- g++ (>= 4.8, for c++11 support)
- swig (>= 3.0)

Note that swig can be installed using the anaconda package manager:

.. code:: bash

   conda install "swig>=3.0"

The CLASS code will automatically be downloaded and compiled, and is thus, not an external dependency for the user. 
The version of CLASS compiled by the code is stored in the variable ``classylss.version.class_version``.

Installation
------------

The package should be compiled using the GNU compilers for C++ and fortran, ``g++`` and ``gfortran``. 
If these are not the default compilers on your system (or if a specific version should be used), they should be 
explicitly set via environment variables. Then, the package can be installed via the `pip` command

.. code:: bash

   # set compilers explicitly, if they are not the default compilers
   export CXX=g++
   export F90=gfortran

   # install the package
   pip install classylss
   
The above procedure has been tested successfully on Mac and Linux machines. However, if installation fails,
it is likely due to a failure while compiling either CLASS or the underlying C++ library. In this case, 
the package should be installed from source and the compilation procedure customized. 

The package can be downloaded from github as

.. code:: bash

    git clone https://github.com/nickhand/classylss.git
    cd classylss

If ``CLASS`` is not built succesfully, the user
can edit the default configuration variables in ``depends/class.cfg``, which are used
when building the ``CLASS`` library.

To verify that the installation has succeeded, run:

.. code-block:: python

    import classylss
    
Examples
--------

To compute power spectra for the Planck 2015 cosmology:

.. code-block:: python

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
    
and similarly, correlation functions can be computed: 

.. code-block:: python
    
    from classylss import correlation
    
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

    
