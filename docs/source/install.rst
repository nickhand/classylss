Installation
============

.. contents::
   :depth: 2
   :local:
   :backlinks: none

classylss is currently supported on macOS and Linux architectures. The
recommended installation method uses
the `Anaconda <https://www.continuum.io/downloads>`_
Python distribution. The package is available for Python versions
2.7, 3.5, and 3.6.


.. _conda-installation:

Installing with Anaconda
------------------------

The easiest installation method uses the ``conda`` utility, as part
of the `Anaconda <https://www.continuum.io/downloads>`_ package
manager. The distribution can be downloaded and installed for free from
https://www.continuum.io/downloads. The package can be installed into the
user's current conda environment using:

.. code:: bash

   conda install -c bccp classylss


Installing from Source
----------------------

.. warning::

    The easiest and recommended method to install classylss
    is using the Anaconda package, which provides pre-built binaries for
    Linux and macOS machines. See :ref:`conda-installation` for more details.

The package can be installed (and compiled locally) via the ``pip`` command

.. code:: bash

   pip install classylss

Users can also downloaded the source code from GitHub using

.. code:: bash

    git clone https://github.com/nickhand/classylss.git
    cd classylss

These installation methods require CLASS to be compiled locally on the user's
machine. If this compilation fails, users can edit the default configuration
variables in ``depends/class.cfg``, which are used when building the CLASS
library. If compilation problems persist, try using the
:ref:`Anaconda installation method <conda-installation>`.

Verifying the Installation
--------------------------

To verify that the installation has succeeded, run:

.. code-block:: python

    import classylss
