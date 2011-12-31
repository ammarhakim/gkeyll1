******************
Decomposed Regions
******************

.. highlight:: c++

Motivation and Overview
-----------------------

To get Lucee to run in parallel an abstraction is needed that stores
the parallel decomposition and provides information needed in
communicating ghost and skin cell values. A class called
``DecompRegion`` is proposed that stores this decomposition and
provides methods for neighbor calculations. A base class providing an
interface for different decomposition algorithms is also proposed.

Proposed Implementation
-----------------------

A *decomposed region* is an abstraction that holds a decomposition of
a given *parent* n-box [#n-box]_ as (potentially) smaller
non-overlapping n-boxes. The union of all the n-boxes in the
decomposition covers the parent box. Such a decomposition is required
for implementing parallel algorithms: the parent box represents the
complete index space of a field while the sub-boxes represent the
portion of the index space handled by individual processors. The class
``Lucee::Region<NDIM, int>`` is used to represent the n-box.::

  template <unsigned NDIM> class DecompRegion 
  {
    public:
  /** Return specified sub-region */
      Lucee::Region<NDIM, int> getRegion(unsigned bn) const;

    private:
  /** Regions making up decomposition */
      std::vector<Lucee::Region<NDIM, int> > rgns;
  };

The decomposed region also allows computing *neighbors* of a given sub
box. The neighbor information, in particular the intersection regions,
can be used for communicating data across processors. The neighbors of
a box are computed by extending it by a specified amount and finding
those boxes that intersect it. Thus, the neighbor calculation can
potentially return different results based on how much the target box
is extended.::

  std::vector<unsigned> getNeigbors(unsigned target, 
    const int lowerExt[NDIM], const int upperExt[NDIM]);

where the arrays ``lowerExt`` and ``upperExt`` specify the number of
cells on the lower and upper faces, respectively, along each
dimension. The aove method returns boxes that share a corner with the
``target`` box. However, a method is also provided that returns only
those boxes that share a common face with the target box::

  std::vector<unsigned> getFaceNeigbors(unsigned target, 
    const int lowerExt[NDIM], const int upperExt[NDIM]);

Different algorithms can be used to compute the decomposition. These
algorithms are implemented a separate set of classes and take the
parent region and the number of sub-boxes as input and compute the
needed decomposition. The base class interface is::

  template <unsigned NDIM> class DecompRegionCalcIfc 
  {
    public:
  /** Compute decomposition, adding subboxes into  decompRgn object */
      virtual void computeDecomp(unsigned nrgns,
        DecompRegion<NDIM>& decompRgn) = 0;
  };

---------------

.. [#n-box] An *n-box* is a rectangular region in n dimensional
   integer lattice space.
