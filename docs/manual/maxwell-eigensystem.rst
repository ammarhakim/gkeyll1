The eigensystem of the Maxwell equations with extension to perfectly hyperbolic Maxwell equations
=================================================================================================

`PDF of note <./_static/files/1012-maxwell-eigsys.pdf>`_

Eigensystem of Maxwell equations
--------------------------------

In this document I list the eigensystem of the Maxwell
equations. Maxwell's equations consist of the curl equations

.. math::

  \frac{\partial \mathbf{B}}{\partial t} + \nabla\times\mathbf{E} &= 0 \\
  \epsilon_0\mu_0\frac{\partial \mathbf{E}}{\partial t} -
  \nabla\times\mathbf{B} &= -\mu_0\mathbf{J}

along with the divergence relations

.. math::

  \nabla\cdot\mathbf{E} &= \frac{\varrho_c}{\epsilon_0} \\
  \nabla\cdot\mathbf{B} &= 0.

Here, :math:`\mathbf{E}` is the electric field, :math:`\mathbf{B}` is
the magnetic flux density, :math:`\epsilon_0`, :math:`\mu_0` are
permittivity and permeability of free space, and :math:`\mathbf{J}`
and :math:`\varrho_c` are specified currents and charges
respectively. The speed of light is determined from
:math:`c=1/(\mu_0\epsilon_0)^{1/2}`.

These are linear equations and hence the eigensytem is independent of
the value of the electromagnetic fields. In 1D Maxwell equations can
be written as, ignoring sources,

.. math::

  \frac{\partial }{\partial t}
  \left[
    \begin{matrix}
      E_x \\
      E_y \\
      E_z \\
      B_x \\
      B_y \\
      B_z
    \end{matrix}
  \right]
  +
  \frac{\partial }{\partial x}
  \left[
    \begin{matrix}
      0 \\
      c^2B_z \\
      -c^2B_y \\
      0 \\
      -E_z \\
      E_y
    \end{matrix}
  \right]
  =
  0.

The eigenvalues of this system are :math:`\{0,0,c,c,-c,-c\}`. The
right eigenvectors of the flux Jacobian are given by the columns of
the matrix

.. math::

  R
  =
  \left[
    \begin{matrix}
      0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 0 & c & 0 & -c & 0 \\
      0 & 0 & 0 & -c & 0 & c \\
      1 & 0 & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & 1 & 0 & 1 \\
      0 & 0 & 1 & 0 & 1 & 0
    \end{matrix}
  \right].

The left eigenvectors are the rows of the matrix

.. math::

  L
  =
  \left[
    \begin{matrix}
      0 & 0 & 0 & 1 & 0 & 0 \\
      1 & 0 & 0 & 0 & 0 & 0 \\
      0 & \frac{1}{2c} & 0 & 0 & 0 & \frac{1}{2} \\
      0 & 0 & -\frac{1}{2c} & 0 & \frac{1}{2} & 0 \\
      0 & -\frac{1}{2c} & 0 & 0 & 0 & \frac{1}{2} \\
      0 & 0 & \frac{1}{2c} & 0 & \frac{1}{2} & 0
    \end{matrix}
  \right].

Eigensystem of Perfectly Hyperbolic Maxwell equations
-----------------------------------------------------

The perfectly hyperbolic Maxwell equations are a modification of the
Maxwell equations that take into account the divergence relations.

.. math::

  \frac{\partial \mathbf{B}}{\partial t} + \nabla\times\mathbf{E} +
  \gamma \nabla\psi
  &= 0 \\
  \epsilon_0\mu_0\frac{\partial \mathbf{E}}{\partial t} -
  \nabla\times\mathbf{B} +     \chi \nabla \phi
  &= -\mu_0\mathbf{J} \\
  \frac{1}{\chi}\frac{\partial \phi}{\partial t} + \nabla\cdot\mathbf{E} 
  &= \frac{\varrho_c}{\epsilon_0} \\
  \frac{\epsilon_0\mu_0}{\gamma}\frac{\partial \psi}{\partial t} + \nabla\cdot\mathbf{B} 
  &= 0.

Here, :math:`\psi` and :math:`\psi` are correction potentials for the
electric and magnetic field respectively and :math:`\chi` and
:math:`\gamma` are dimensionless factors that control the speed at
which the errors are propagated.

In 1D these equations can be written as, ignoring sources,

.. math::

  \frac{\partial }{\partial t}
  \left[
    \begin{matrix}
      E_x \\
      E_y \\
      E_z \\
      B_x \\
      B_y \\
      B_z \\
      \phi \\
      \psi
    \end{matrix}
  \right]
  +
  \frac{\partial}{\partial x}
  \left[
    \begin{matrix}
      \chi c^2 \phi \\
      c^2B_z \\
      -c^2B_y \\
      \gamma \psi \\
      -E_z \\
      E_y \\
      \chi E_x \\
      \gamma c^2B_x
    \end{matrix}
  \right]
  =
  0.

The eigenvalues of this system are :math:`\{-c\gamma, c\gamma, -c\chi,
c\chi, c, c, -c, -c\}`. The right eigenvectors of the flux Jacobian
are given by the columns of the matrix

.. math::
  R
  =
  \left[
    \begin{matrix}
      0  & 0 & 1 & 1 & 0 &  0 &  0 & 0 \\
      0  & 0 & 0 & 0 & c &  0 & -c & 0 \\
      0  & 0 & 0 & 0 & 0 & -c &  0 & c \\
      1  & 1 & 0 & 0 & 0 &  0 &  0 & 0 \\
      0  & 0 & 0 & 0 & 0 &  1 &  0 & 1 \\
      0  & 0 & 0 & 0 & 1 &  0 &  1 & 0 \\
      0  & 0 & -\frac{1}{c} & \frac{1}{c} & 0 &  0 &  0 & 0 \\
     -c  & c & 0 & 0 & 0 &  0 &  0 & 0
    \end{matrix}
  \right].

The left eigenvectors are the rows of the matrix

.. math::
  L
  =
  \left[
    \begin{matrix}
      0 & 0 & 0 & \frac{1}{2} & 0 & 0 & 0 & -\frac{1}{2c} \\
      0 & 0 & 0 & \frac{1}{2} & 0 & 0 & 0 & \frac{1}{2c} \\
      \frac{1}{2} & 0 & 0 & 0 & 0 & 0 & -\frac{c}{2} & 0 \\
      \frac{1}{2} & 0 & 0 & 0 & 0 & 0 & \frac{c}{2} & 0 \\
      0 & \frac{1}{2c} & 0 & 0 & 0 & \frac{1}{2} & 0 & 0 \\
      0 & 0 & -\frac{1}{2c} & 0 & \frac{1}{2} & 0 & 0 & 0 \\
      0 & -\frac{1}{2c} & 0 & 0 & 0 & \frac{1}{2} & 0 & 0 \\
      0 & 0 & \frac{1}{2c} & 0 & \frac{1}{2} & 0 & 0 & 0
    \end{matrix}
  \right].
