==================================================
Lucee: Light-weight Universal Computational Engine
==================================================

Lucee programs are written in the `Lua <http://www.lua.org>`_
programming language, version 5.2. Most valid Lua programs are also
valid Lucee programs [#lua-except]_. The Lua programming language is
described in the book `Programming in Lua, Second Edition
<http://www.inf.puc-rio.br/~roberto/pil2>`_. You need not buy this
book: the first edition of the book is available for free `here
<http://www.lua.org/pil>`_. The firt edition is sufficient for almost
all Lucee simulations.

In this manual the Lucee specific Lua objects and methods are
described. A Lucee simulation is created by writing a Lua program that
uses these Lucee-specific objects along with standard Lua control
structures and functions. This gives great flexibility as a powerful
general-purpose language is available to you to create a simulation.

Contents
--------

.. toctree::
  :maxdepth: 2

  datastruct
  hyperequation

.. [#lua-except] Importing external modules is not presently
   supported. The reason for this is that most supercomputers do not
   allow loading shared libraries, a required feature for external Lua
   modules to work.
