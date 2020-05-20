****************************
Laminate Theory Calculator
****************************

Introduction
============
This set of scrips contains functions that can be used for the design of simple
fibre composite laminates. The goal of this project is to create a open-source
laminate theory calculator something that can compete with "The Laminator".
There are two benefits of the project in python is that it can easily be
integrated in, for example, a stiffend plate buckling/cripling calculator.
Or that it can be used to find the optimal layout by iterating over various
designs.

Before using the code I recommend to read upon literature introducing the
Classical Laminate Theory. Many such books exist, from my experience the
following where reccomended:

- Design and Analysis of Composite Structures: With Applications to Aerospace Structures by Christos Kassapoglou (ISBN: 9781118536933, DOI: `10.1002/9781118536933 <https://doi.org/10.1002/9781118536933>`_)
- Composite Materials: Design and Applications by Daniel Gay (ISBN: 9780429101038 DOI: `10.1201/b17106 <https://doi.org/10.1201/b17106>`_)

.. todo::

   - Adding functions to the homonigization module.
   - Adding a inverse ply failure calculator for Tsai-Hill.
   - Adding buckling / cripling calculators.

.. toctree::
   :caption: Table of Contents
   :maxdepth: 1

   homogenization
   abdcal
   deformation
   failure
   license


Example
=======
The example disscussed below is formed by snipplets of `example.py <https://github.com/AJJLagerweij/Classical-Composite-Laminate-Theory-Calculator/blob/master/example.py>`_.

At the start of the file the required packages are needed. A minimum requirement
is `numpy` and the different scripts related to this project

.. literalinclude:: /../../example.py
   :language: python
   :lines: 13-19
   :lineno-match:

Currently one is required to define the ply parameters in the ply axis system.
In the future this will be replaced by a script where a basic homonigization is
performed on the constituants of each ply.

.. literalinclude:: /../../example.py
   :language: python
   :lines: 26-43
   :lineno-match:

Then the laminate properties musth be defined. Starting with the stacking
sequence which consists of a list of the rotation angles of each ply (global to
ply axis system) and a list with the thickness of each ply and the ply stiffness
matrix :math:`Q`. All list must be orderd from the top to the bottom of the
laminate. Notice that the positive :math:`z` direction is downward by convention.

Afterwards the ply properties can be used to calculate the ABD matrix and its
inverse.

.. literalinclude:: /../../example.py
   :language: python
   :lines: 49-56
   :lineno-match:

Now a load or deformation vector can be applied. Here the load vector was used.
The load vector is a 1 by 3 numpy matrix consists of the running loads and
moments in the form of :math:`(N_x, N_y, N_{xy}, M_x, M_y, M_{xy})^T`.
Similarly a deformation vector is defined as
:math:`(\varepsilon_x, \varepsilon_y, \varepsilon_{xy},\kappa_x, \kappa_y, \kappa_{xy})^T`.
Afterwards the resulting ply stresses and loads (in their local axis system) 
can be calculated. The strain and stress calculations are performed at the top
and bottom of each ply. Detials can be found in the documentation of the
deformation module.

.. literalinclude:: /../../example.py
   :language: python
   :lines: 62-68
   :lineno-match:

Lastly the stresses are used to calculate if the failure criterias are
violated. Here the the max stress, Tsai-Wu and Tsai-Hill criterias are used.
It is reccomended that the user reads up on the differences between the possible
criteria, all of them have their specific strength and weaknesses and are meant
for their specif purpose. If one does not keep this in mind properly one will
end up with flawed designs.

.. literalinclude:: /../../example.py
   :language: python
   :lines: 74-77
   :lineno-match:

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

