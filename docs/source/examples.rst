Examples
========

.. contents::
    :depth: 3
    :local:
    :backlinks: none

.. ipython:: python
    :suppress:

    import matplotlib.pyplot as plt
    plt.style.use('notebook.mplstyle')

.. note::

    For a full description of the available functionality, see the :ref:`API`
    section.


The main binding of CLASS is available to users in the :mod:`classylss.binding`
module. To start, we shall make the necessary imports and print out the
:mod:`classylss` version and CLASS version:

.. ipython:: python

    import classylss
    import classylss.binding as CLASS

    print("classylss version = ", classylss.__version__)
    print("CLASS version = ", classylss.class_version)

Initializing the CLASS engine
-----------------------------

Users must first initialize the CLASS engine, optionally passing in a dictionary
of parameters. If no parameters are passed, the default CLASS parameters are used.
This can be done using the :class:`classylss.binding.ClassEngine` class as

.. ipython:: python

    engine = CLASS.ClassEngine({'H0':70, 'Omega_m':0.31})
    default_engine = CLASS.ClassEngine() # default parameters

Loading parameters from file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can use the :func:`classylss.load_ini` function to load parameters
from a CLASS ``.ini`` file into a dictionary. For example, we can load
the `concise.ini <https://cdn.rawgit.com/lesgourg/class_public/master/concise.ini>`_
parameter file from the CLASS GitHub and use these parameters to initialize
our :class:`~classylss.binding.ClassEngine`:

.. ipython:: python

  params = classylss.load_ini('concise.ini')
  print(params)

  engine = CLASS.ClassEngine(params)

Loading precision parameters from file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can also pass in precision parameters that have been loaded from file
using the :func:`classylss.load_precision` function. For example, we can load
the `pk_ref.pre <https://cdn.rawgit.com/lesgourg/class_public/master/pk_ref.pre>`_
parameter file from the CLASS GitHub:

.. ipython:: python

  pre_params = classylss.load_precision('pk_ref.pre')
  print(pre_params)

  # default cosmo params + precision params
  high_pre_engine = CLASS.ClassEngine(pre_params)


Using the main CLASS modules
----------------------------

The Background module
~~~~~~~~~~~~~~~~~~~~~

The `background module <https://goo.gl/SU71dn>`_ in CLASS computes
background cosmology quantities, e.g., distances, as a function of redshift.
It also provides access to the parameters of the specified cosmological model.
In :class:`classylss`, the :class:`classylss.binding.Background` object provides
an interface to this module.

.. ipython:: python

    # initialize the background module
    bg = CLASS.Background(engine)

    # print out some cosmological parameters
    print("h = ", bg.h)
    print("Omega0_m = ", bg.Omega0_m)
    print("Omega0_lambda = ", bg.Omega0_lambda)
    print("Omega0_r = ", bg.Omega0_r)
    print("Omega0_k = ", bg.Omega0_k)

Distances
^^^^^^^^^

We can plot distances measures as a function of redshift using:

.. ipython:: python

  import numpy as np
  from matplotlib import pyplot as plt

  z = np.linspace(0., 2., 1024)

  plt.plot(z, bg.angular_diameter_distance(z), label=r"$D_A$")
  plt.plot(z, bg.comoving_distance(z), label=r"$D_C$")
  plt.plot(z, bg.luminosity_distance(z), label=r"$D_L$")

  # save
  plt.legend()
  plt.xlabel(r"$z$")
  @savefig distances.png
  plt.ylabel(r"distance $[h^{-1} \mathrm{Mpc}]$")


Linear growth evolution
^^^^^^^^^^^^^^^^^^^^^^^

The linear growth factor and growth rate can be computed as:

.. ipython:: python

  z = np.linspace(0., 2., 1024)

  plt.plot(z, bg.scale_independent_growth_factor(z), label=r"$D(z)$")
  plt.plot(z, bg.scale_independent_growth_rate(z), label=r"$f(z)$")

  # save
  plt.legend()
  plt.xlabel(r"$z$")
  @savefig growth.png
  plt.ylabel("growth factor/rate")

Density parameter evolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The density parameters :math:`\Omega` can be computed as a function of
redshift for various species:

.. ipython:: python

  z = np.linspace(0., 2., 1024)

  plt.loglog(1+z, bg.Omega_m(z), label=r"$\Omega_m$")
  plt.loglog(1+z, bg.Omega_r(z), label=r"$\Omega_r$")

  # save
  plt.legend()
  plt.xlabel(r"$1+z$")
  @savefig density-parameters.png
  plt.ylabel("density parameter")

Setting ``a_max``
^^^^^^^^^^^^^^^^^

Quantities in the background module can also be computed for :math:`a > 1` by
setting the ``a_max`` parameter. Below, we compute the same distance calculations
up to :math:`a_\mathrm{max}=2`:

.. ipython:: python

  a_max = 2.0
  cosmo = CLASS.ClassEngine({'a_max':a_max})
  bg = CLASS.Background(cosmo)

  z = np.linspace(1/a_max-1, 1, 1024)

  plt.plot(z, bg.angular_diameter_distance(z), label=r"$D_A$")
  plt.plot(z, bg.comoving_distance(z), label=r"$D_C$")
  plt.plot(z, bg.luminosity_distance(z), label=r"$D_L$")

  # save
  plt.legend()
  plt.xlabel(r"$z$")
  @savefig a-max-distances.png
  plt.ylabel(r"distance $[h^{-1} \mathrm{Mpc}]$")


The Spectra module
~~~~~~~~~~~~~~~~~~

The `spectra module <https://goo.gl/EMti1s>`_ in CLASS computes
linear density and velocity transfer functions, as well as the linear
and nonlinear density power spectra. In :class:`classylss`, the
:class:`classylss.binding.Spectra` object provides an interface to this module.

.. ipython:: python

  # initialize with proper output
  cosmo = CLASS.ClassEngine({'output': 'dTk vTk mPk', 'non linear': 'halofit', 'P_k_max_h/Mpc' : 20., "z_max_pk" : 100.0})
  sp = CLASS.Spectra(cosmo)


The linear and nonlinear power spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can compute the nonlinear and linear power spectra at the desired redshifts
and wavenumbers using:

.. ipython:: python

  z = 0.5
  k = np.logspace(-2, 0, 100)

  # nonlinear power
  pk_nl = sp.get_pk(k=k, z=z)

  # linear power
  pk_lin = sp.get_pklin(k=k, z=z)

  plt.loglog(k, pk_nl, label='nonlinear')
  plt.loglog(k, pk_lin, label='linear')

  # save
  plt.legend()
  plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
  @savefig linear-nonlinear-power.png
  plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^3]$")

The evolution of :math:`\sigma_8`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The evolution of the perturbation amplitude :math:`\sigma_8` can be computed
using:

.. ipython:: python

  z = np.linspace(0., 2., 1024)
  plt.plot(1 + z, sp.sigma8_z(z))

  # save
  plt.xlabel(r"$1+z$")
  @savefig sigma8_z.png
  plt.ylabel(r"$\sigma_8(z)$")

Density and velocity transfers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can compute the transfer functions in CLASS units for density, velocity and
:math:`\phi` in CLASS units, using :func:`~classylss.binding.Spectra.get_transfer`:

.. ipython:: python

    # transfer at z=0
    transfer = sp.get_transfer(z=0)
    print(transfer.dtype.names)

    plt.subplot(211)
    plt.plot(transfer['k'], transfer['d_tot'])
    plt.ylabel("total density transfer")

    plt.subplot(212)
    plt.plot(transfer['k'], transfer['t_tot'])
    plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
    @savefig transfers.png
    plt.ylabel("total velocity transfer")

The Thermo module
~~~~~~~~~~~~~~~~~

The `thermo module <https://goo.gl/JKGUP6>`_ in CLASS computes
various quanitites related to the thermodynamic history of the Universe.
In :class:`classylss`, the :class:`classylss.binding.Thermo` object provides
access to several of the computed quantites, including:


.. ipython:: python

    th = CLASS.Thermo(cosmo)

    print("drag redshift = ", th.z_drag)
    print("sound horizon at z_drag = ", th.rs_drag)
    print("reionization optical depth = ", th.tau_reio)
    print("reionization redshift = ", th.z_reio)
    print("recombination redshift = ", th.z_rec)
    print("sound horizon at recombination = ", th.rs_rec)
    print("sound horizon angle at recombination = ", th.theta_s)

The Primordial module
~~~~~~~~~~~~~~~~~~~~~

The `primordial module <https://goo.gl/SmxLQz>`_ in CLASS computes
various quanitites related to the initial, primordial conditions of the Universe.
In :class:`classylss`, the :class:`classylss.binding.Primordial` object provides
access to the primordial power spectrum, which is defined as:

.. math::

    \mathcal{P_R}(k) = A_s \left (\frac{k}{k_0} \right )^{n_s - 1 + 0.5 \ln(k/k_0) (dn_s / d\ln k)}.


Below, we compute this quantity from the :class:`~classylss.binding.Primordial`
class, and compare to the above equation:

.. ipython:: python

    sp = CLASS.Spectra(cosmo)
    pm = CLASS.Primordial(cosmo)

    k = np.logspace(-3, 0, 100)
    plt.loglog(k, pm.get_pk(k), label='from CLASS')
    plt.loglog(k, sp.A_s * (k / sp.k_pivot)**(sp.n_s-1), ls='--', c='k', label='analytic')

    plt.legend()
    plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
    @savefig primordial-power.png
    plt.ylabel("dimensionless power")
