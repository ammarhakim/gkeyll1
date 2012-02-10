*************************************
Datastructures: Module ``DataStruct``
*************************************

Datastructures are fundamental object required to perform
simulations. Most datastructures live on grid objects and allow the
storage of multiple values per location in the grid.

.. contents::

Fields: ``DataStruct.Field1D``, ``DataStruct.Field2D`` and ``DataStruct.Field3D``
=================================================================================

The blocks ``DataStruct.Field1D``, ``DataStruct.Field2D`` and
``DataStruct.Field3D`` can be used to create 1D, 2D and 3D fields
[#seven-field]_. Fields are fundamental datastructures that allow
attaching data to a grid. A field element corresponds to a location in
the grid and allows the storage of ``numComponents`` number of
elements at that location.

Constructor Parameters
----------------------

In the following table parameters for only 3D table constructor are
listed. The parameters for 1D and 2D fields are identical.

.. list-table:: ``DataStruct.Field3D``
  :header-rows: 1
  :widths: 30,10,60

  * - Variable [Units]
    - Default
    - Description
  * - onGrid
    - None
    - Name of grid object on which this field lives
  * - numComponents
    - 1
    - Number of elements to store at each field location
  * - ghost
    - {0, 0}
    - Number of ghost cells along each side of field.
  * - location
    - "center"
    - Location of data. Should be one of "vertex" or "center".

**Notes** The parameter ``ghost`` is a 2 element table (irrespective
of field dimension) that indicates the number of ghost cells to
create. For example the table ``{1, 2}`` will create 1 ghost cell on
the *lower* edge in each direction and 2 ghost cells on the *upper*
edge in each direction.

The parameter ``location`` specifies the location inside the cell at
which data is assumed to be stored. The data can be stored at the cell
centroid or the lower vertex. This is used only in initializing the
field.

Methods
-------

The 1D, 2D and 3D field object support the following methods.

.. py:function:: write(flNm)

  Write out data in the field to file named *flNm*. The file is
  created if it does not exist and overwritten if it does. Depending
  on the file extension, the data can be written as either HDF5 or
  plain-text [#plain-txt]_.

.. py:function:: clear(val)

  Set all field values to *val*.

.. py:function:: copy(fld)

  Copy from supplied field, *fld*. The field to copy from must live on
  the same grid and have the same number of components as the target
  field.

.. py:function:: duplicate()

  Creates and returns a duplicate of this field.

.. py:function:: alias(s, e)

  Creates and returns an alias of this field. The alias field points
  to :math:`[s,e)` components of the original field. I.e, the first
  component of the alias refers to the s-th component of the original
  field. No new memory is allocated: the alias and the original share
  memory, and hence modifying one will change the other.

.. py:function:: accumulate(coeff_1, fld_1, coeff_2, fld_2, ...)

  Add :math:`\sum coeff_i*fld_i` to the values in this field. Here
  *coeff_i* is a number and *fld_i* is a field that lives on the same
  grid and has the same number of components as the target field.

.. py:function:: combine(coeff_1, fld_1, coeff_2, fld_2, ...)

  Set field to :math:`\sum coeff_i*fld_i`. Here *coeff_i* is a number
  and *fld_i* is a field that lives on the same grid and has the same
  number of components as the target field.

.. py:function:: set(luaFunc)

  This method takes a Lua function to initialize the field. The
  function *luaFunc* must take in the :math:`(x,y,z)` coordinates
  (irrespective of field dimension) and return ``numComponents``
  values, one for each component of the field.

Examples
--------

.. code-block:: lua

  eulerEqn = HyperEquation.Euler {
   gasGamma = 1.4,
  }

.. [#seven-field] Lucee supports the creation of upto seven
   dimensional fields. However, not all implemented algorithms work
   with fields with dimensionality higher than three.

.. [#plain-txt] Writing to plain text is not a good idea except for
   quick debugging. The writes are very slow and do not work in
   parallel.
