.. BSG documentation master file, created by
   sphinx-quickstart on Mon Jan 22 17:51:32 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The BSG and NME Decay Libraries
===============================

The BSG library is a C++ program designed to calculate the allowed beta spectrum shape coming from nuclear beta decay to high precision. It couples to a smaller library, NME, which calculates nuclear matrix elements required for the higher order corrections to the spectrum shape. A graphical user interface is provided to automate parts of the procedure.

The formalism is described in this publication_, while the code is published here_.

.. _publication: https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.90.015008
.. _here: https://www.sciencedirect.com/science/article/pii/S0010465519300645

Installation & Basic execution
------------

.. toctree::
   :maxdepth: 2

   installation

Overview
--------

.. toctree::
   :maxdepth: 2

   intro
   input_files
   spectrum_mod
   nme_options

API Reference
-------------

.. toctree::
   :maxdepth: 2

   API_reference

Status
------

.. include:: status.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
