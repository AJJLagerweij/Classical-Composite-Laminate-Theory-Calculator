**********************
Deformation Calculator
**********************

Introduction
============
This file containts functions to calculate what the deformation or load on a
laminate is. To perform this calculation the ABD matrix and its inverse are
required in combination with either:

- Load vector :math:`(N_x, N_y, N_{xy}, M_x, M_y, M_{xy})^T`
- Deformation vector :math:`(\varepsilon_x, \varepsilon_y, \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`

Eventually it will calculate the stress and strain values inside each ply of the
laminate.

.. contents:: Table of Contents
   :depth: 2

Routines
========
.. automodule:: deformation
   :members:

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

