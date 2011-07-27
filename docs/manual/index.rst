==================================================
Lucee: Light-weight Universal Computational Engine
==================================================

Lucee programs are written in the `Lua <http://www.lua.org>`_
programming language, version 5.2. Most valid Lua programs are also
valid Lucee programs [#lua-except]_. The Lua programming language is
described in the book `Programming in Lua, Second Edition
<http://www.inf.puc-rio.br/~roberto/pil2>`_. You need not buy this
book: the first edition of the book is available for free `here
<http://www.lua.org/pil>`_. The first edition is sufficient for almost
all Lucee simulations.

In this manual the Lucee specific Lua objects and methods are
described. A Lucee simulation is created by writing a Lua program that
uses these Lucee-specific objects along with standard Lua control
structures and functions. This gives great flexibility as a powerful
general-purpose language is available for creating a simulation.

Lua accessible object documentation
-----------------------------------

.. toctree::
  :maxdepth: 2

  mathphys
  datastruct
  hyperequation

Technical notes
---------------

These technical notes give details, mainly mathematical, about
algorithms and physical models implemented in Lucee. Various parts of
the code refer to these notes which should be considered as references
for the implemented equations and schemes.

.. toctree::
  :maxdepth: 2

  maxwell-eigensystem
  euler-eigensystem
  hancock-muscl

.. [#lua-except] Importing external modules is not presently
   supported. The reason for this is that most supercomputers do not
   allow loading shared libraries, a required feature for external Lua
   modules to work. Also, allowing use of arbitrary libraries reduces
   reproducibility. If you really want a module you will need to add
   to in the source code and recompile Lucee.
