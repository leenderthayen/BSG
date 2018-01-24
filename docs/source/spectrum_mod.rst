.. _spectrum_mod_label:

Using the BSG spectrum calculation
==================================

The beta spectrum of the outgoing beta particle or (anti)neutrino consists of a great number of terms due to the complex environment it is created in. An in-depth discussion on the physical origin and mathematical derivations and discussions can be found in the original work (arXiv_). A summary and additional background on the workings of this code can be found here_.

.. _arXiv: https://arxiv.org/abs/1709.07530
.. _here: https://arxiv.org

Each correction term can be turned on and off individually. The following table shows all included corrections with their corresponding notation and expected magnitude.

.. image:: table_spectral_corrections.png

Toggling corrections
--------------------

Each correction term can be turned on and off individually. The user has two options to do this:

- Using the command line
- Using the configuration file, discussed in :ref:`input_files_label`.

An overview of these options can be found by running the program's ``--help`` command, where among other things one can find

.. literalinclude:: spectral_options_help.txt

The simple options are to be used on the command line, while those between brackets can be added to the configuration file using the typical INI style, with the ``Spectrum`` header.

Advanced options
----------------

The options described in the previous section should be plentiful for most users. In special cases, however, it might be necessary to introduce further customisation into some of the corrections. All of these are related to the explicit shape of the (weak) charge distribution of the nucleus, and as such are related to two things

- The electrostatic higher-order finite size correction, *U*
- The convolution finite size correction, *C*

Electrostatic finite size correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When moving from the point-charge nucleus to a uniformly charged sphere, we have added the L0 correction the the typical F0 Fermi function. As those is often too simplistic a picture, the *U* correction was introduced. In the original work (arXiv_), two possibilities were given which are implemented already in the code, chosen using the ``Spectrum.ESShape`` option 

- Fermi distribution, chosen using the ``Fermi`` option
- Modified Gaussian distribution, chosen using the ``Modified Gaussian`` option

If the user requires a custom potential, the code allows for a specification of the the first 3 terms in a power expansion of the electrostatic potential. As an example, the potential for a uniformly charged sphere is given by

.. math::

   V(r) = -\frac{\alpha Z}{R}(\frac{3}{2}-\frac{r^2}{2R^2})

We define now the coefficients :math:`v_n` so that

.. math::
   V(r) = \sum_{n=0}^\infty\left(-\frac{\alpha Z}{R^{n+1}} \right)v_nr^n

so that for the uniformly charged sphere we have :math:`v_0=3/2`, :math:`v_2=-1/2` and :math:`v_4=0`. If the user wants to define a new *U* correction using a different potential, one can set these options for the new potential as such. The old :math:`v_n` are then those of the potential against one wishes to correct. In the typical case, this would be those of the uniformly charged sphere. 

   Note that due to the limited number of terms in the original formulation, this new correction is limited in precision to lower *Z* nuclei.

Convolution finite size correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just like the user has the option to choose the electrostatic charge distribution, so can one do the same for the convolution finite size correction. In a first approximation, the *C* correction is calculated assuming the weak charge density and the simple charge density distributions to be one and the same. For the latter then, the user has the same options as for the electrostatic counterpart in the previous section, specified through the ``Spectrum.CShape`` correction.

As the weak charge is typically not the same as the simple charge distribution, an additional correction was defined: the isovector correction, :math:`C_I`.
