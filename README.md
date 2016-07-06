# Gkeyll: Computational Plasma Physics Framework

Gkeyll programs are written in the Lua (http://www.lua.org)
programming language, version 5.2. Most valid Lua programs are also
valid Gkeyll programs. The Lua programming language is
described in the book "Programming in Lua, Second Edition"
http://www.inf.puc-rio.br/~roberto/pil2. You need not buy this
book: the first edition of the book is available for free at
http://www.lua.org/pil. The first edition is sufficient for almost
all Gkeyll simulations.

A Gkeyll simulation is created by writing a Lua program that
uses these Gkeyll-specific objects along with standard Lua control
structures and functions. This gives great flexibility as a powerful
general-purpose language is available for creating simulations.

For extensive documentation of various problems performed with Gkeyll, please see http://ammar-hakim.org/sj/.

You need to create a Bitbucket account and then you can clone the repo on you local machine. We use the Mercurial version control system, and you need to download and install that before you clone Gkeyll.

If you want to modify Gkeyll code, you will need permission set in Bitbucket. Only the project admins can do that. Please contact either Ammar Hakim, Bhuvana Srinivasan or John Loverich if you think you need write access to the code.


# Getting and Installing Gkeyll

If you are working at PPPL, University of New Hampshire, Virginia Tech or have an account on NERSC or the University of Texas Stampede cluster, you do not need to build Gkeyll. In fact, you should not build your own copy. Gkeyll is installed at these locations and you can use the pre-built versions.

If you must build Gkeyll, please remember that the process may be complicated, although in principle the whole system can be built with a *single command*.