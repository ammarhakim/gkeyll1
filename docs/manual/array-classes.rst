Array classes
-------------

Lucee provides several array classes to store data in indexed
containers. The array classes can be instantiated for any type,
although some algorithms and functions are defined only for ``float``
and ``double``.

``Lucee::FixedVector``: Fixed-size vectors
++++++++++++++++++++++++++++++++++++++++++

The class ``Lucee::FixedVector`` provides a one-dimensional,
fixed-size array. The size of the array must be known at compile
time. For example, to create a vector with three double-precision
elements one can use::

  Lucee::FixedVector<3, double> v1(1.0);
  Lucee::FixedVector<3, double> v2(1.0, 2.0, 3.0);

This will create ``v1`` and initialize all its elements to 1.0 and
create ``v2`` and initialize its elements to 1.0, 2.0 and 3.0. The
elements can be accessed using standard array notation::

  double sum = 0.0;
  for (unsigned i=0; i<3; ++i)
    sum += v1[i]*v1[i];

For the complete API of this class see the Doxygen generated
documentation for `Lucee::FixedVector`_.

.. _Lucee::FixedVector: ../../../api/html/class_lucee_1_1_fixed_vector.html

``Lucee::Array``: :math:`N`-dimensional arrays
++++++++++++++++++++++++++++++++++++++++++++++

.. class:: Array

  The class ``Lucee::Array`` provides a :math:`N`-dimensional
  reference counted array to store data of arbitrary types. It is
  declared as::

    namespace Lucee
    {
      template <unsigned N, typename T> class Array;
    }

  The first template parameter indicates the rank of the array and the
  second, the type of data stored in the array.

  .. cfunction:: ctor Array (unsigned shape[N], const T& init)
    :noindex:

    Create a new array with given shape. Start indices are assumed to
    be :math:`(0,\ldots)`. An optional ``init`` value can be specified
    and is applied to all elements of the array. By default
    ``init=0``.

  .. cfunction:: ctor Array (unsigned shape[N], int start[N], const T& init)
    :noindex:

    Create a new array with given shape and specified start
    indices. An optional ``init`` value can be specified and is
    applied to all elements of the array. By default ``init=0``.

  .. cfunction:: ctor Array (const Array<T>& arr)
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

  .. cfunction:: bool isContiguous ()
    :noindex:

    Return true if the array is stored contiguously. Trying to create
    views or accessing the raw memory pointer of a non-contiguous
    array will lead to a run-time exception being thrown.

  .. cfunction:: void fillWithShape (unsigned shape[N])
    :noindex:

    On return fill the shape of the array in ``shape``.

  .. cfunction:: void fillWithStart (int start[N])
    :noindex:

    On return fill the start index of the array in ``start``.

  .. cfunction:: unsigned getShape (unsigned dir)
    :noindex:

    Return the shape of the array in direction ``dir``.

  .. cfunction:: int getLower (unsigned dir)
    :noindex:

    Return the starting index in direction ``dir``.

  .. cfunction:: int getUpper (unsigned dir)
    :noindex:

    Return *one past* the last index in direction ``dir``.

  .. cfunction:: T& operator() (int i, int j, ...)
    :noindex:

    Access element at :math:`(i,j,\ldots)` index location in
    array. For a :math:`N`-dimensional array exactly N indices
    must be specified. A compile-time error will occur when trying to
    use this method with the incorrect number of indices.

  .. cfunction:: T& first ()
    :noindex:

    Return reference to the first element in array. This is useful
    when passing the raw pointer to the array data to functions
    expecting ``T*``. An exception is thrown if the array is not
    contiguous.

  .. cfunction:: T& operator() (int i[N])
    :noindex:

    Access element at index :math:`(i_1,\ldots,i_N)` index location in
    array.

  .. cfunction:: Array<N,T> createView (unsigned shape[N] , int start[N], int newStart[N])
    :noindex:

    Returns a view into the original array. The portion of the
    original array indexed by the view-array is specified by lower
    bounds, ``start``, and the shape, ``shape``. The indices
    ``newStart`` indicate the new starting index of the view-array. In
    most cases this can simply be set to ``start``.

``Lucee::Vector``: One-dimensional arrays
+++++++++++++++++++++++++++++++++++++++++

.. class:: Vector

  The ``Lucee::Vector`` class inherits from the ``Lucee::Array`` class
  and provides a specialized container for one-dimensional arrays of
  ``int``, ``float`` or ``double``. Certain linear algebra functions
  are only defined for ``float`` or ``double``. It is declared as::

    namespace Lucee
    {
      template <typename T> class Vector : public Array<1, T>;
    }

  .. cfunction:: ctor Vector (unsigned len)
    :noindex:

    Create a new vector with specified length. The start index is
    assumed :math:`0`. All vector elements are initialized to :math:`0`.

  .. cfunction:: ctor Vector (unsigned len, int start)
    :noindex:

    Create a new vector with specified length and start index
    ``start``. All vector elements are initialized to :math:`0`.

  .. cfunction:: T& operator[] (int i)
    :noindex:

    Access element at index :math:`i` in vector.

  .. cfunction:: unsigned getLength ()
    :noindex:

    Get length of the vector.

  .. cfunction:: Vector<T> duplicate ()
    :noindex:

    Duplicate the vector. The returned vector is contiguous and has
    identical data as this vector.

``Lucee::Matrix``: Matrix class
+++++++++++++++++++++++++++++++

.. class:: Matrix

  The ``Lucee::Matrix`` class inherits from the ``Lucee::Array`` class
  and provides a specialized container for two-dimensional arrays of
  ``int``, ``float`` or ``double``. Certain linear algebra functions
  are only defined for ``float`` or ``double``. The matrix elements
  are stored in column major order to enable use of Fortran routines
  for linear algebra. It is declared as::

    namespace Lucee
    {
      template <typename T> class Matrix : public Array<2, T, Lucee::ColMajorIndexer<2> >
    }

  .. cfunction:: ctor Matrix (unsigned row, unsigned col)
    :noindex:

    Create a new matrix with specified ``row`` and ``col``. The start
    indices are assumed :math:`(0,0)`. All matrix elements are
    initialized to :math:`0`.

  .. cfunction:: ctor Matrix (unsigned shape[2])
    :noindex:

    Create a new matrix with specified ``shape``. The start indices
    are assumed :math:`(0,0)`. All matrix elements are initialized to
    :math:`0`.

  .. cfunction:: ctor Matrix (unsigned shape[2], int start[2])
    :noindex:

    Create a new matrix with specified ``shape`` and given ``start``
    indices. All matrix elements are initialized to :math:`0`.

  .. cfunction:: Matrix<T> duplicate ()
    :noindex:

    Duplicate the matrix. The returned matrix is contiguous and has
    identical data as the original matrix.

  .. cfunction:: bool isSquare ()
    :noindex:

    Is the matrix square?

  .. cfunction:: bool isTranspose ()
    :noindex:

    Is the matrix transpose of another matrix?

  .. cfunction:: unsigned numRows ()
    :noindex:

    Return number of rows in matrix.

  .. cfunction:: unsigned numColumns ()
    :noindex:

    Return number of columns in matrix.

  .. cfunction:: Vector<T>& getColumn (unsigned col)
    :noindex:

    Return the ``col`` column of the matrix.

  .. cfunction:: Vector<T>& getRow (unsigned row)
    :noindex:

    Return the ``row`` row of the matrix.

  .. cfunction:: Matrix<T> transpose ()
    :noindex:

    Create the transpose of the matrix. No data is actually allocated
    and the transposed matrix shares data with the original
    matrix. Hence, modifying one will affect the other.

Linear algebra functions
++++++++++++++++++++++++

Linear algebra functions are defined in the header file
``LcLinAlgebra.h``. The available linear algebra methods are described
below. In the following the template type is either ``double`` or
``float``.

.. cfunction:: Matrix<double>& accumulate (double beta, Matrix<double>& C, double alpha, const Matrix<double>& A, const Matrix<double>& B)

  Compute the matrix-matrix product :math:`C \leftarrow \alpha AB +
  \beta C`. An exception is thrown if the :math:`A` and :math:`B`
  matrices are not of the correct shape. A reference to :math:`C` is
  returned.

.. cfunction:: Matrix<double>& accumulate (Matrix<double>& C, const Matrix<double>& A, const Lucee::Matrix<double>& B)

  Compute the matrix-matrix product :math:`C \leftarrow AB`.  An
  exception is thrown if the :math:`A` and :math:`B` matrices are not
  of the correct shape. A reference to :math`C` is returned.

.. cfunction:: Vector<double>& accumulate (double beta, Vector<double>& y, double alpha, const Matrix<double>& A, const Vector<double>& x)

  Compute the matrix-vector product :math:`y \leftarrow \alpha Ax +
  \beta y`. An exception is thrown if the :math:`A`, :math:`x` and
  :math:`y` are not of the correct shape. A reference to :math:`y` is
  returned.

.. cfunction:: Vector<double>& accumulate (Vector<double>& y, const Matrix<double>& A, const Vector<double>& x)

  Compute the matrix-vector product :math:`y \leftarrow Ax`. An
  exception is thrown if the :math:`A`, :math:`x` and :math:`y` are
  not of the correct shape. A reference to :math:`y` is returned.

.. cfunction:: Matrix<double>& accumulate (Matrix<double>& A, double alpha, const Vector<double>& x, const Vector<double>& y)

  Compute the vector-vector outer product :math:`A \leftarrow \alpha
  xy^T + A`. An exception is thrown if the :math:`A`, :math:`x` and
  :math:`y` are not of the correct shape. A reference to :math:`A` is
  returned.

.. cfunction:: void eig (const Matrix<T>& A, Vector<T>& evr, Vector<T>& evi)

  Compute the eigenvalues of the matrix :math:`A`. The matrix must be
  square or an exception is thrown. The real part of the eigenvalues
  are returned in ``evr`` and the imaginary part are returned in
  ``evi``. These vectors must be pre-allocated and contiguous.

.. cfunction:: void eig (const Matrix<T>& A, Vector<T>& evr, Vector<T>& evi, Matrix<T>& vecl, Matrix<T>& vecr)

  Compute the eigenvalues and the eigenvectors of the matrix
  :math:`A`. The matrix must be square or an exception is thrown. The
  real part of the eigenvalues are returned in ``evr`` and the
  imaginary part are returned in ``evi``. The left eigenvectors are
  returned as columns of ``vecl`` and the right eigenvectors as the
  columns of ``vecr``. The vectors and matrices must be pre-allocated
  and contiguous.

.. cfunction:: void eigRight (const Matrix<T>& A, Vector<T>& evr, Vector<T>& evi, Matrix<T>& vec)

  Compute the eigenvalues and the right eigenvectors of the matrix
  :math:`A`. The matrix must be square or an exception is thrown. The
  real part of the eigenvalues are returned in ``evr`` and the
  imaginary part are returned in ``evi``. The right eigenvectors are
  returned as the columns of ``vec``, which must be of the same shape
  as this matrix. The vectors and matrices must be pre-allocated and
  contiguous.

.. cfunction:: void eigLeft (const Matrix<T>& A, Vector<T>& evr, Vector<T>& evi, Matrix<T>& vec)

  Compute the eigenvalues and the left eigenvectors of the matrix. The
  matrix must be square or an exception is thrown. The real part of
  the eigenvalues are returned in ``evr`` and the imaginary part are
  returned in ``evi``. The left eigenvectors are returned as the
  columns of ``vec``, which must be of the same shape as this
  matrix. The vectors and matrices must be pre-allocated and
  contiguous.

.. cfunction:: void solve (const Matrix<T>& A, Matrix<T>& rhs)

  Solve the linear system of equations :math:`Ax=b`, where :math:`b`
  are columns of the matrix ``rhs``. The matrix ``A`` must be square
  or an exception is thrown. The ``rhs`` matrix must have same number
  of rows as ``A``. On output the columns of the ``rhs`` matrix are
  replaced by the corresponding solution vector. A LU-decomposition is
  used to solve the system of equations.

``Lucee::ColMajorIndexer``: Column major indexer
++++++++++++++++++++++++++++++++++++++++++++++++

.. class:: ColMajorIndexer

  The ``Lucee::ColMajorIndexer`` class provides a mapping of a
  N-dimensional index space into a 1-dimensional linear
  index space. It is declared as::

    namespace Lucee
    {
      template <unsigned NDIM> class ColMajorIndexer;
    }

  The template parameter ``NDIM`` indicates the rank (dimension) of
  the index space.

  .. cfunction:: ctor ColMajorIndexer (unsigned shape[NDIM], int start[NDIM])
    :noindex:

    Create a indexer which maps a ``NDIM`` dimensional index space of
    specified ``shape`` and given ``start`` indices into a linear
    1-dimensional index.

  .. cfunction:: int getLower (unsigned dir)
    :noindex:

    Return the starting index in direction ``dir``.

  .. cfunction:: int getUpper (unsigned dir)
    :noindex:

    Return *one past* the last index in direction ``dir``.

  .. cfunction:: int getIndex (int i, int j, ...)
    :noindex:

    Return index into the 1-dimensional space corresponding to the
    index in the N-dimensional space :math:`(i,j,\ldots)`. For a
    N-dimensional space exactly N indices must be specified. A
    compile-time error will occur when trying to use this method with
    the incorrect number of indices.

  .. cfunction:: int getGenIndex (int i[N])
    :noindex:

    Return index into the 1-dimensional space corresponding to the
    index :math:`(i_1,\ldots,i_N)` index location in N-dimensional
    space.

