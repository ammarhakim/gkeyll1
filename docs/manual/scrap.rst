Array classes
-------------

To start a new input file start off with::

  tStart = 0.0
  tEnd = 1.0
  numFrames = 10
  dtInit = 1.0e-4

  <warpx>
    kind = fluidComponent

    <grid>
      kind = cartGrid1D
      lower = [0, 0]
    </grid>
  </warpx>

A common use of arrays in Lucee is to store vector fields. For example,
to store the electric field in a 2D simulation one would create a 3D
array (the final index for the components of the electric field)::

  unsigned shape[3] = {16, 32, 3};
  Lucee::Array<3, double> E(shape);

Although this array can be indexed in the usual way, it is ofter more
convenient to get access to all three components::

  Lucee::ArrayItr<double> ef = E.createItr();
  E.setItr(ef, 3, 4);
  ef[0] = 1.0; ef[1] = 2.0; ef[2] = 1.0;

The iterator ``ef`` points to the three elements of the electric
field at location :math:`(3,4)`. This code is equivalent to::

  E(3, 4, 0) = 1.0;
  E(3, 4, 1) = 2.0;
  E(3, 4, 2) = 1.0;

The advantage of using the iterator to access the components is that the
iterator can be passed to other functions which expect ``double *``::

 double norm = computeNorm(&ef[0]);

where, for example,::

  double computeNorm(double *vec) 
  {
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  }


Maxwell's Equations
+++++++++++++++++++

Of course, one must remember that the Maxwell equations are

.. math::

   \pfrac{f}{t} + v \frac{\partial f}{\partial x} = 0

It that cool or what? Now, one can take the moments of this equation
to get

.. math::
  :nowrap: 

  \notag
  \begin{align}
    a &= b \\
    c &= d
  \end{align}

So AMSMATH should work with. Below are Maxwell equations.

.. math::
  :nowrap:

  \notag
  \begin{align}
    \nabla\times \mvec{E} &= -\pfrac{\mvec{B}}{t} \\
    \nabla\times \mvec{B} &=
    \mu_0\mvec{J}+\frac{1}{c^2}\pfrac{\mvec{E}}{t} \\
    \nabla\cdot\mvec{E} &= \frac{\varrho_c}{\varepsilon_0} \\
    \nabla\cdot\mvec{B} &= 0.
  \end{align}

Did this work or not? So here is one more test.

``Lucee::Array``: Basic array class
+++++++++++++++++++++++++++++++++++

.. class:: Array

  This is the base array class. It is declared as::

    template <typename T> class Array;

.. cfunction:: int Array<T>. Array(unsigned shape[NDIM], int start[NDIM])
   :noindex:

   A ctor.

.. cfunction:: int Array<T>. Array(unsigned shape[NDIM])
   :noindex:

   Another ctor.

.. cfunction:: int Array<T>. Array(unsigned shape[NDIM], int start[NDIM])
   :noindex:

   Yet another ctor.

.. cfunction:: unsigned Array<T>. getDims(int x)
   :noindex:

   This is a C function.

.. cfunction:: T& Array. getData(int x)
   :noindex:

   This is another C function.

.. cfunction:: int makeArray(int *s, int *e)

   Create a new array.