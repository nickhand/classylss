classylss
----------

.. image:: https://img.shields.io/pypi/v/classylss.svg
   :alt: PyPi
   :target: https://pypi.python.org/pypi/classylss/

.. image:: https://api.travis-ci.org/nickhand/classylss.svg
    :alt: Build Status
    :target: https://travis-ci.org/nickhand/classylss/

A lightweight Python binding of the CMB Boltzmann code `CLASS`, with an emphasis on the routines that are important for large-scale structure calculations. The main modules of the CLASS code are exposed to the user via a Cython wrapper.

.. _`CLASS` : http://class-code.net


Dependencies
------------

The package is lightweight and the only dependencies are:

- numpy
- cython

The CLASS code will automatically be downloaded and compiled, and is thus, not an external dependency for the user. However, the user will need a valid C compiler to compile the CLASS code. The version of CLASS compiled by the code is stored in the variable ``classylss.class_version``.

Installation
------------

The package can be installed as part of the Anaconda package manager using

.. code:: bash

   conda install -c bccp classylss

The package can also be installed via the `pip` command

.. code:: bash

   pip install classylss
   
The package can be also be downloaded from github using

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

See the tests of the code in ``classylss/tests/`` for examples of using each of the main CLASS modules. 
   
Feature Requests
----------------

Additional features of the CLASS code that are not yet implemented can be exposed via the Cython wrapper relatively easily. We encourage users to open up a GitHub issue with any feature requests.
