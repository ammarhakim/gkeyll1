Array classes
-------------

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

