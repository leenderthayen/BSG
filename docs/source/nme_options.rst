.. _nme_options_label:

Nuclear Matrix Element usage
============================

The ``NME`` library calculates nuclear matrix elements in the Behrens-Bühring formalism, meaning these are noted as follows

.. math::
   ^{V/A}\mathcal{M}_{KLs}^{(n)}

Currently only the first order in :math:`n` is allowed, i.e., :math:`n=0` and does not need to be specified.

The command line options are as follows:

- ``-b``: Calculate the weak magnetism contribution normalized with the Gamow-Teller form factor and mass number according to Holstein, i.e., :math:`b/Ac_1`.
- ``-d``: Calculate the first-class induced tensor contribution, likewise normalized with the Gamow-Teller form factor and mass number according to Holstein, i.e., :math:`d/Ac_1`
- ``-M V/AKLs``: Calculate the general matrix element :math:`^{V/A}\mathcal{M}_{KLs}^{(0)}`

Currently, only matrix elements appearing in allowed :math:`\beta` decay are supported

The nuclear matrix element due to some operator :math:`\mathcal{O}` is most easily trated when writing it in second quantization

.. math::

   M_i = \langle f | \mathcal{O} | i \rangle = \sum_{\alpha \beta} \langle \alpha | \mathcal{O} | \beta \rangle \langle f | a^\dagger_\alpha a_\beta | i \rangle

Here :math:`\alpha, \beta` are single particle states, and :math:`\langle f | a^\dagger_\alpha a_\beta | f \rangle` are the one-body density matrix elements (OBDME), calculated in nuclear many-body calculations such as the shell model.

As many matrix elements or form factors depend on the value of the axial-vector coupling constant, :math:`g_A`, yet it is often quenched to values below 1.2723, the configuration file introduces a new constant ``Constants.gAeff`` that can be set to its proper value.

The code gives two possibilities for the calculation of these matrix elements:

- The simple approach based on the extreme single-particle approximation
- The more advanced method of loading in OBDME from many-body calculations

Extreme single-particle
-----------------------

In this case, we consider only a single initial and final single state, meaning the sum in the previous equation disappears. All responsability now lies with the correct designation of these states. The program offers several options for this endeavor. All of these are specified under the ``Computational`` header in the general configuration file, defaulted to ``config.txt``.

- ``Method``: The method by which to calculate the nuclear matrix elements. By default this will be ``ESP``: Extreme single-particle.
- ``Potential``: Sets the nuclear potential. We have three possibilities:

   + ``SHO``: The spherical harmonic oscillator potential. In this case, nuclear single-particle states are pure harmonic oscillator functions as determined by filling nucleons in the regular *jj*-coupling.
   + ``WS``: The spherical Woods-Saxon potential, with spin-orbit coupling. Single-particle wave functions are combinations of spherical harmonic oscillator functions and typically correspond nicely to the *j*-coupling results.
   + ``DWS``: The deformed Woods-Saxon potential, with spin-orbit coupling. As *j* is now any more a good quantum number, we consider the projection along the axial symmetry axis, *K*. Wave functions are combinations of spherical harmonic oscillator functions for all :math:`j \geq K`.

The potential for the (deformed) Woods-Saxon potential is as follows

.. math::

   \mathcal{H} = -\frac{\hbar^2}{2m}\nabla^2 - V_0f(r) - V_s\left(\frac{\hbar}{m_\pi c}\right)^2\frac{1}{r}\frac{df}{dr}l\cdot s + V_0R\frac{df}{dr}\sum^6_{n\text{ even}}(\beta_{n}Y^0_n)

where :math:`V_0` is the depth of the Woods-Saxon potential, :math:`f(r) = (1+\exp[(r-R)/a_0])^{-1}` is the Woods-Saxon function with surface thickness :math:`a_0` and nuclear radius :math:`R`, and spherical harmonics, :math:`Y_{L}^M`. The spin-orbit term contains the coupling constant, :math:`V_s`, and the Compton wavelength of the pion, :math:`\hbar/m_\pi c`. In the case of protons, an additional Coulomb potential is added. The depth of the Woods-Saxon potential is given in its 'optimized' form as

.. math::
   V_0 = V\left(1\pm \chi \frac{N-Z}{N+Z}  \right)

with the upper (lower) sign for protons (neutrons), and :math:`V=-49.6` MeV and :math:`\chi=0.86` by default.

These parameters can be changed separately for protons and neutrons according to the following options, present in the configuration file

- ``Vproton/neutron``: The initial depth of the Woods-Saxon potential
- ``Xproton/Xneutron``: The asymmetry factor :math:`\chi` in the potential depth
- ``V0Sproton/neutron``: The depth of the spin-orbit potential
- ``SurfaceThickness``: The surface thickness, :math:`a_0` in femtometer.

As the simple filling schemes based on these potentials do not always accurately predict the correct valence spin state, the code allows for an enforcement of the correct single-particle states within a user-specified energy margin. It contains several options to enforce spin selection and coupling

- ``ForceSpin``: Instead of choosing the single-particle state that is obtained through simple filling, pick the closest one corresponding to the spin state defined in the ``Mother/Daughter.ForcedSPSpin`` in even-*A* nuclei. In odd-*A* nuclei, the correct spin state is automatically chosen if one can be found within the user-defined energy window.
- ``ReversedGallagher``: Most examples of deformed, even-*A* nuclei follow the simple spin selection rules by Gallagher. In the original work, a 'reversed' selection rule was defined, which can be turned on in the code.
- ``OverrideSPCoupling``: In case a properly coupled state cannot be obtained using conventional coupling rules, override the coupling completely and set the coupled spins to the corresponding nuclear state.
- ``EnergyMargin``: Set the margin in MeV to select a different state corresponding to the proper initial or final state when it is not obtained through regular methods.

When activated, the state closest in energy to the originally proposed level with the correct spin state is selected as the one that participates in the interaction.

Once the state is selected, matrix elements are calculated and output is written to a ``.nme`` file, containing information on the wave function composition of initial and final states and the calculation results. An example excerpt is given below

.. literalinclude:: excerpt_output_31S.nme

Many-body input
---------------

While the extreme single-particle can achieve good results with minimal effort, given a good basis, it is oftentimes necessary to use results from nuclear many-body calculations such as the nuclear shell model. Currently, there is support for this by specifying the single-particle state :math:`\alpha` and :math:`\beta` defined above together with their OBMDE in a CSV file as follows

.. code-block:: bash

   j_i, n_i, l_i, j_f, n_f, l_f, OBDME

in a spherical harmonic oscillator basis.
