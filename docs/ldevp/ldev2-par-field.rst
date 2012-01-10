**********************
LDEV2: Parallel Fields
**********************

.. highlight:: c++

.. contents::

Motivation and Overview
-----------------------

To run algorithms in parallel we need a distributed data-structure
that allows storing different portions of a field on different
processors. To enable this an algorithm is required that splits the
domain based on the number of processors. Further, a method is needed
that synchronizes data shared between processes. The first of these
requirements is already met by the ``DecompRegion`` class and the
algorithms that create objects of this class. A ``sync()`` method will
be added to the ``Field`` class that performs the synchronization.

Proposed Implementation
-----------------------

Lua use case
++++++++++++

Before discussing the proposed implementation it is important to see
how parallel fields are created and used from Lua. Here is an
example

.. code-block:: lua

  -- create decomposition object
  decomp = DecompRegionCalc2D.CartProd { cuts = {2, 2} }

  grid = Grid.RectCart2D {
     lower = {-0.5, -0.5},
     upper = {0.5, 0.5},
     cells = {50, 50},
     decomposition = decomp, -- decomposer to use
  }

  q = DataStruct.Field2D {
     onGrid = grid,
     -- [Bx, By, Ez, psi]
     numComponents = 4,
     ghost = {1, 1},
     decompose = true, -- this 'true' by default
  }

In this code snippet we first create the object ``decomp`` that
performs the decomposition. In this case we are using the cartesian
product decomposer that allows specifying the number of "cuts" to use
in each direction. The grid is uses the key ``decomposition`` to
indicate which decomposer the grid should use to distribute itself
over multiple processors. Multiple grids can use the same
decomposition object. A field created over this grid uses this
decomposition to create the needed parallel data-structure.

Implementation
++++++++++++++

To implement the needed code two classes need to be modified: the
``StructGridBase`` class needs to hold the decomposition object while
the ``StructGridField`` class needs to implement the ``sync()``
method. Essentially, the ``StructGridField`` object gets a hold of the
decomposition object from the grid it lives on and uses the neighbor
information it provides to send/receive the data it needs for
communication. The ``sync()`` method needs to be a pure virtual in the
base class ``Field`` to allow for future field types (multi-block
fields, unstructured fields, etc.) to implement their specific
synchronization methods.
