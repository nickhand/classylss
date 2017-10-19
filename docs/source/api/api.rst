.. _api:

API Reference
=============

The CLASS engine can be initialized from a dictionary of parameters using:

.. autosummary::

    ~classylss.binding.ClassEngine


We provide wrappers for the main modules in the CLASS code:

.. autosummary::

    ~classylss.binding.Background
    ~classylss.binding.Perturbs
    ~classylss.binding.Primordial
    ~classylss.binding.Spectra
    ~classylss.binding.Thermo

There is also a wrapper object that provides compatibility with the syntax
used by the :mod:`astropy` cosmology objects:

.. autosummary::

    ~classylss.astropy_compat.AstropyCompat

And we provide two functions for loading parameters from file into a dictionary
(which can be passed to :class:`~classylss.binding.ClassEngine`):

.. autosummary::

    classylss.load_ini
    classylss.load_precision
