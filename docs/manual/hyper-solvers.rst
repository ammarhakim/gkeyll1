Algorithms for solution of hyperbolic balance laws
--------------------------------------------------

An important class of algorithms implemented in Lucee are for the
solution of multidimension hyperbolic balance laws. These laws can be
written as

.. math::
  :nowrap:

  \notag
  \begin{align}
    \pfrac{q}{t} + \pfrac{f}{x} + \pfrac{g}{y} + \pfrac{h}{z} = s
  \end{align}

where :math:`q(x,y,z,t)` is a system of m conserved variables,
:math:`f(q), g(q), h(q)` are flux functions and :math:`s(q,t)` are
sources. This set of equations is called a system of *hyperbolic
balance laws* if the Jacobian matrix

.. math::
  :nowrap:

  \notag
  \begin{align}
    A(q,n) = n_x \pfrac{f}{q} + n_y \pfrac{g}{q} + n_z \pfrac{h}{q}
  \end{align}

where :math:`n = (n_x,n_y,n_z)` are components of a unit vector, has
real eigenvalues and a complete set of right eigenvectors. Further,
the system is called *isotropic* if the eigenvalues do not depend on
:math:`(n_x,n_y,n_z)`. In general, the solution to such systems
preserve *invariant domains*, i.e. solutions satisfy :math:`q \in
\script{U}`, where :math:`\script{U}` is an open subset of
:math:`R^m`. For example, for neutral fluid equations the energy and
density must be non-negative.

There are three major classes of solvers for hyperbolic balance laws
in Lucee: the *wave-propagation scheme*, the *MUSCL scheme* and the
*nodal discontinuous Galerkin* scheme. Each of these schemes works in
general geometries and preserves the invariant domains for that
system.