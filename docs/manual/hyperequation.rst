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
:math:`\varepsilon(p,\rho) = p/\rho(\gamma-1)`.

The object assumes that the conserved variables are stored in the
order :math:`[\rho, \rho u, \rho v, \rho w, E]`, where :math:`\rho` is
the fluid density, :math:`(u,v,w)` are the components of the fluid
velocity and :math:`E` is the fluid total energy given by

.. math::

  E = \rho \varepsilon + \frac{1}{2}\rho (u^2+v^2+w^2)

Note that all three components of the momentum (or velocity) are
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
  * - minPressure [Pa]
    - 0.0
    - Minimum allowable pressure (see text for details)
  * - minDensity [Kg/m :math:`^3`]
    - 0.0
    - Minimum allowable density (see text for details)

**Notes** The optional variables ``minDensity`` and ``minPressure``
represent the minimum pressure and density computed by the Euler
equation object. These parameters are useful to prevent the formation
of negative density or pressure.

Methods
-------

This object does not support any methods.

Examples
--------

.. code-block:: lua

  eulerEqn = HyperEquation.Euler {
   gasGamma = 1.4,
  }

Implementation notes
--------------------

The ``HyperEquation.Euler`` object represents the Euler equation
written in the conservative form. It pro **give link to doxygen**

