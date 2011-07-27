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

The ``HyperEquation.Euler`` block creates an ideal Euler equation
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

The following code block uses the euler equation block and conserved
variables field ``qCons`` to compute the primitive variables ``qPrim``
and then write out the pressure to an HDF5 file

.. code-block:: lua

  eulerEqn = HyperEquation.Euler { gasGamma = 1.4 }

  -- compute primitive variables
  eulerEqn:primitive(qCons, qPrim)

  -- create alias to point to pressure
  pressure = qPrim:alias(4, 5) -- pressure is 4th component (density is 0th)

  -- write out pressure
  pressure:write("pressure.h5")

---------

Maxwell equations: ``HyperEquation.Maxwell``
============================================

The ``HyperEquation.Maxwell`` block creates an Maxwell equation object
that represents the Maxwell equations of electromagnetism in free
space.

.. math::

  \frac{\partial \mathbf{B}}{\partial t} + \nabla\times\mathbf{E} &= 0 \\
  \epsilon_0\mu_0\frac{\partial \mathbf{E}}{\partial t} -
  \nabla\times\mathbf{B} &= 0

Here, :math:`\mathbf{E}` is the electric field, :math:`\mathbf{B}` is
the magnetic flux density, :math:`\epsilon_0`, :math:`\mu_0` are
permittivity and permeability of free space.

The object assumes that the fields are stored in the order
:math:`[E_x, E_y, E_z, B_x, B_y, B_z]`. Note that all field components
are stored even for 1D and 2D simulations.

The object created with this block can be used in updaters that solve
hyperbolic equations.

Constructor Parameters
----------------------

.. list-table:: ``HyperEquation.Maxwell``
  :header-rows: 1
  :widths: 30,10,60

  * - Variable [Units]
    - Default
    - Description
  * - lightSpeed [m/s]
    - None
    - Speed of light in free space.

Methods
-------

The Maxwell equation object does not support and special methods.

Examples
--------

To create an Maxwell equation object with speed of light taken from
Lucee defined values

.. code-block:: lua

  maxwellEqn = HyperEquation.Maxwell {
   lightSpeed = Lucee.SpeedOfLight,
  }

---------

Perfectly Hyperbolic Maxwell equations: ``HyperEquation.PhMaxwell``
===================================================================

The ``HyperEquation.PhMaxwell`` block creates a Perfectly Hyperbolic
Maxwell equation object that represents the Maxwell equations of
electromagnetism in free space.

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

Here, :math:`\mathbf{E}` is the electric field, :math:`\mathbf{B}` is
the magnetic flux density, :math:`\epsilon_0`, :math:`\mu_0` are
permittivity and permeability of free space. Also, :math:`\psi` and
:math:`\psi` are correction potentials for the electric and magnetic
field respectively and :math:`\chi` and :math:`\gamma` are
dimensionless factors that control the speed at which the errors are
propagated.

The object assumes that the fields are stored in the order
:math:`[E_x, E_y, E_z, B_x, B_y, B_z, \phi, \psi]`. Note that all
field components are stored even for 1D and 2D simulations.

The object created with this block can be used in updaters that solve
hyperbolic equations.

Constructor Parameters
----------------------

.. list-table:: ``HyperEquation.Maxwell``
  :header-rows: 1
  :widths: 30,10,60

  * - Variable [Units]
    - Default
    - Description
  * - lightSpeed [m/s]
    - None
    - Speed of light in free space.
  * - elcErrorSpeedFactor
    - 0.0
    - Value of :math:`\chi`. Error propagation speed is :math:`\chi c`
  * - mgnErrorSpeedFactor
    - 0.0
    - Value of :math:`\gamma`. Error propagation speed is :math:`\gamma c`

Methods
-------

The Maxwell equation object does not support and special methods.

Examples
--------

To create an Perfectly Hyperbolic Maxwell equation object with speed
of light taken from Lucee defined values

.. code-block:: lua

  maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = Lucee.SpeedOfLight,
   elcErrorSpeedFactor = 1.0,
   mgnErrorSpeedFactor = 1.0
  }