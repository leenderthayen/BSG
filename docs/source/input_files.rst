.. _input_files_label:

Input files
===========

Both the ``BSG`` and ``NME`` library function using two input files

- Transition-specific input file
- General configuration file

Both are written in typical Windows-INI style. While both libraries add a number of options to both of these, there is significant overlap. We discuss these here in turn. Options specific to both libraries can be found at :ref:`spectrum_mod_label` and :ref:`nme_options_label`.

Transition-specific input file
------------------------------

This file deals with the beta transition for which the quantities are to be calculated. It contains three headers: ``Transition``, ``Mother`` and ``Daughter``. The last two are identical in form. An example is given below for the :math:`\beta^+` decay of :math:`^{19}` Ne

.. literalinclude:: 19Ne.ini
   :language: ini

General configuration file
--------------------------

The general configuration file defaults to ``config.txt`` and contains straightforward options. A documented example is given below

.. literalinclude:: config_overlap.txt
   :language: ini
