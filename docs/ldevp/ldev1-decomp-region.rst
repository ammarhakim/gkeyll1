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

Decomposed Region
+++++++++++++++++

A *decomposed region* is an abstraction that holds a decomposition of
a given *global* region [#region]_ as (potentially) smaller
non-overlapping regions. The union of all the regions in the
decomposition covers the global region. Such a decomposition is
required for implementing parallel algorithms: the global region
represents the complete index space of a field while the sub-regions
represent the portion of the index space handled by individual
processors. The class ``Lucee::Region<NDIM, int>`` is used to
represent the region.::

  template <unsigned NDIM> class DecompRegion 
  {
    public:
  /** Return specified sub-region */
      Lucee::Region<NDIM, int> getRegion(unsigned bn) const;

    private:
  /** Regions making up decomposition */
      std::vector<Lucee::Region<NDIM, int> > rgns;
  };

The decomposed region also allows computing *neighbors* of a given
sub-region. The neighbor information, in particular the intersections,
can be used for communicating data across processors. The neighbors of
a region are computed by extending it by a specified amount and
finding those regions that intersect it. Thus, the neighbor
calculation can potentially return different results based on how much
the target region is extended.::

  std::vector<unsigned> getNeighbors(unsigned target, 
    const int lowerExt[NDIM], const int upperExt[NDIM]);

where the arrays ``lowerExt`` and ``upperExt`` specify the number of
cells on the lower and upper faces, respectively, along each
dimension. The above method returns, in addition to face neighbors,
regions that share a corner with the ``target`` region. However, a
method is also provided that returns only those regions that share a
common face with the target region::

  std::vector<unsigned> getFaceNeighbors(unsigned target, 
    const int lowerExt[NDIM], const int upperExt[NDIM]);

Depending on the stencil of the algorithm the complete neighbor
information (including corners) might not be needed for communication
of ghost and skin data.

Computing the Decomposition
+++++++++++++++++++++++++++

Different algorithms can be used to compute the decomposition. These
algorithms are implemented a separate set of classes and take the
global region and the number of sub-regions as input and compute the
needed decomposition. The base class for the decomposition algorithm
supports the interface::

  template <unsigned NDIM> class DecompRegionCalcIfc 
  {
    public:
  /** Calculate decomposition, adding subregions into object */
      void calcDecomp(unsigned nrgns, DecompRegion<NDIM>& decompRgn);
  };

This method is passes an existing decomposition with the global region
already set. It decomposes the region into the specified number of
sub-regions, modifying the ``decompRgn`` object with this information.

Derived classes that implement the decomposition algorithm must
provide the method::

  virtual void decompose(unsigned nrgns, 
    const Region<NDIM, int>& globalRgn) = 0;

Here, the ``globalRgn`` is the global region that must be decomposed
into ``nrgns`` sub-regions. Once it computes the sub-regions this
method should add each of them by successive calls to the base class
method::

  void addRegion(const Region<NDIM, int>& subRgn);

The base class ensures that the current decomposition in
``decompRegion`` is cleared out and that the decomposition provided by
the derived classes covers ``globalRgn``.

---------------

.. [#region] An *region* is a rectangular region in n dimensional
   integer lattice space.
