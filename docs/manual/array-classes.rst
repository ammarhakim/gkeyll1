Array classes
-------------

Lucee provides several array classes to store data efficiently.

``Lucee::FixedVector``: Fixed-size vectors
++++++++++++++++++++++++++++++++++++++++++

.. class:: Lucee::FixedVector

  The class ``Lucee::FixedVector`` provides a one-dimensional array of
  fixed size. The size of the array must be know at compile time. It
  is declared as::

    namespace Lucee
    {
      template <unsigned NELEM, typename T>
      class FixedVector;
    }

  The first template parameter specifies the size of the vector and
  second parameter specifies the type of data stored in the vector.

  .. cfunction:: ctor FixedVector (const T& init)
    :noindex:

    Create a new fixed-size vector and set all values to ``init``.

  .. cfunction:: ctor FixedVector (const T& init)
    :noindex:

    Create a new fixed-size vector and set all values to ``init``.

  .. cfunction:: ctor FixedVector (T vals[NELEM])
    :noindex:

    Create a new fixed-size vector and set values to ones specified in
    the ``val`` array.

  .. cfunction:: ctor FixedVector (T v1, ...)
    :noindex:

    Create a new fixed-size vector and set values :math:`(v_1,
    \ldots)`. Note that exactly ``NELEM`` values must be specified for
    this function to work correctly. The method will fail silently if
    incorrect number of values are specified.

  .. cfunction:: T& operator[](int i)
    :noindex:

    Access element at index :math:`i` in vector.

  .. cfunction:: unsigned getLength ()
    :noindex:

    Get number of elements in vector.

  .. cfunction:: T getNorm ()
    :noindex:

    Compute :math:`l_2`-norm of vector. The norm only makes sense for
    double and floating point vectors.

``Lucee::Array``: :math:`N`-dimensional arrays
++++++++++++++++++++++++++++++++++++++++++++++

.. class:: Lucee::Array

  The class ``Lucee::Array`` provides a :math:`N`-dimensional
  reference counted array to store data of arbitrary types. It is
  declared as::

    namespace Lucee
    {
      template <unsigned NDIM, typename T> class Array;
    }

  The first template parameter indicates the rank of the array and the
  second, the type of data stored in the array.

  .. cfunction:: ctor Array (unsigned shape[NDIM], const T& init)
    :noindex:

    Create a new array with given shape. Start indices are assumed to
    be :math:`(0,\ldots)`. An optional ``init`` value can be specified
    and is applied to all elements of the array. By default
    ``init=0``.

  .. cfunction:: ctor Array (unsigned shape[NDIM], int start[NDIM], const T& init)
    :noindex:

    Create a new array with given shape and specified start
    indices. An optional ``init11 value can be specified and is
    applied to all elements of the array. By default ``init=0``.

  .. cfunction:: ctor Array (const Array& arr)
    :noindex:

    Create a new array that is a shallow copy of ``arr``. No data is
    actually allocated and the new array points to the same memory
    space as ``arr``.

  .. cfunction:: Lucee::Array& operator= (const Array<T>& arr)
    :noindex:

    Shallow copy ``arr``. The call creates an alias for ``arr``,
    i.e. no data is allocated and modifying one changes the other.

  .. cfunction:: Lucee::Array& operator= (const T& val)
    :noindex:

    Set all elements of array to ``val``.

  .. cfunction:: unsigned getRank ()
    :noindex:

    Return the rank of the array.

Array indexing and sequencing
+++++++++++++++++++++++++++++

When an array is created the data is stored in a single chunk of
contiguous memory. Hence, a method of mapping of the
:math:`N`-dimensional index space to a :math:`1`-dimensional index
space is needed. There are two mapping functions provided in Lucee:
row-major and column-major indexing. These are implemented in the
``Lucee::ColMajorIndexer`` and the ``Lucee::RowMajorIndexer`` classes.

Let :math:`(i_1,\ldots,i_N)` be the index in the :math:`N`-dimensional
index space. Let :math`(s_0,\ldots,s_N)` be the starting index and
:math:`(l_0,\ldots,l_N)` be the shape of the space. Then, a linear
mapping, :math:`0\le k<L`, where :math:`L=\Pi_{n=0}^N l_n`, can be
constructed as

.. math::

  k(i_1,\ldots,i_N) = a_0 + \sum_{n=1}^N a_n i_n

where :math:`a_n`, :math:`i=1,\ldots,N` are coefficients determined by the
particular indexing method, and

.. math::

  a_0 = -\sum_{n=1}^N a_n s_n.

In column-major indexing the region in the :math:`1`-dimensional space
spanned by the last index is contigous:

.. math::

  k(i_1,\ldots,s_N+1) = k(i_1,\ldots,s_N) + 1

Using this condition in the mapping function yields the recursion
relation :math:`a_{j-1}=a_j l_j` with :math:`a_N=1`.

In row-major indexing the region in the :math:`1`-dimensional space
spanned by the first index is contigous:

.. math::

  k(s_1+1,\ldots,i_N) = k(s_1,\ldots,i_N) + 1

Using this condition in the mapping function yields the recursion
relation :math:`a_{j+1}=a_j l_j` with :math:`a_1=1`.
