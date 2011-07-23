**********************************************
Hyperbolic Equations: Module ``HyperEquation``
**********************************************

Hyperbolic equation objects are defined in the ``HyperEquation``
module. The constructors in this module can be used to create equation
objects that can be used in various solvers that evolve hyperbolic
equations. These objects also allow calculating primitive variables
from conserved variables.

.. contents::

Euler equations: ``HyperEquation.Euler``
========================================

The ``HyperEquation.Euler`` block create an ideal Euler equation
object that represents a gas with a constant adiabatic index
:math:`\gamma`. I.e, the internal energy of the gas is computed from
:math:`\varepsilon(p,\rho) = p/\rho(\gamma-1)`, where :math:`p` is the
gas pressure and :math:`\rho` is the density.

The object assumes that the conserved variables are stored in the
order :math:`[\rho, \rho u, \rho v, \rho w, E]`. The
``HyperEquation.Euler`` object decribes the 3D Euler equation written
in the form

.. math::
  :label: eq:euler-eqn

  \frac{\partial}{\partial{t}}
  \left[
    \begin{matrix}
      \rho \\
      \rho u \\
      \rho v \\
      \rho w \\
      E
    \end{matrix}
  \right]
  +
  \frac{\partial}{\partial{x}}
  \left[
    \begin{matrix}
      \rho u \\
      \rho u^2 + p \\
      \rho uv \\
      \rho uw \\
      (E+p)u
    \end{matrix}
  \right]
  =
  0

Here only the X-direction fluxes are show. Here, :math:`(u,v,w)` are
the components of the fluid velocity and :math:`E` is the fluid total
energy given by

.. math::

  E = \rho \varepsilon + \frac{1}{2}\rho (u^2+v^2+w^2)

Note that all three components of the momentum (and velocity) are
stored even for 1D and 2D simulations.

The object created with this block can be used in
updaters that solve hyperbolic equations.

Constructor Parameters
----------------------

.. list-table:: ``HyperEquation.Euler``
  :header-rows: 1
  :widths: 30,10,60

  * - Variable [Units]
    - Default
    - Description
  * - gasGamma
    - None
    - Gas adiabatic index :math:`\gamma`

Methods
-------

The Euler equation object supports the following methods.

.. py:function:: primitive(cons, prim)

  Given conserved variables stored in the *cons* field, computes and
  stores the primitive variables in the *prim* field. The primitive
  variables are stored in the order :math:`[\rho, u, v, w, p]`, where
  :math:`p` is the fluid pressure.

.. py:function:: conserved(prim, cons)

  Given primitive variables stored in the *prim* field, computes and
  stores the conserved variables in the *cons* field.

Examples
--------

To create an Euler equation object with :math:`\gamma = 1.4` you can do

.. code-block:: lua

  eulerEqn = HyperEquation.Euler {
   gasGamma = 1.4,
  }
