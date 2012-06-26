SciMake CONVENTIONS FOR CMAKE FILES AND MODULES

$Id: README.txt 1162 2011-12-17 16:42:43Z cary $

The txcmake document is at
https://ice.txcorp.com/support/wiki/CmakeCodingStandards.
This needs to be rewritten for the new structure.

Links for conventions
http://www.cmake.org/pipermail/cmake/2004-July/005283.html
http://www.phy.bnl.gov/~bviren/lbne/code/ai/external/build/LCG/cmake-2.6.4/Modules/readme.txt

To work and commit in this directory, you will need to

  svn switch --relocate svn://svn.code.sf.net/p/scimake/code/trunk https://SOURCEFORGE_USERNAME@svn.code.sf.net/p/scimake/code/trunk

# Older solutions:
  mv scimake scimake.sav
  svn co --username=SOURCEFORGE_USERNAME svn+ssh://SOURCEFORGE_USERNAME@svn.code.sf.net/p/scimake/code/trunk scimake

If we change the external, one could do
  svn co http://svn.code.sf.net/p/scimake/code/trunk scimake
  cd scimake
  svn switch --relocate http://svn.code.sf.net/p/scimake/code/trunk https://SOURCEFORGE_USERNAME@svn.code.sf.net/p/scimake/code/trunk



