Array indexing and sequencing
-----------------------------

This note describes the indexing and sequencing algorithms implemented
in Lucee.

Array indexing
++++++++++++++

An *indexing function* maps an :math:`N`-dimensional index space to a
:math:`1`-dimensional index space. Indexing functions are used in the
Lucee array classes to map array indices into a linear index. There
are several sets of mapping classes provided in Lucee: row-major,
column-major and completely symmetric-tensor indexing. These are
implemented in the ``Lucee::ColMajorIndexer`` and the
``Lucee::RowMajorIndexer`` classes. 

Let :math:`(i_1,\ldots,i_N)` be the index in the :math:`N`-dimensional
index space. Let :math:`(s_1,\ldots,s_N)` be the starting index and
:math:`(l_1,\ldots,l_N)` be the shape of the space. Then, a linear
mapping, :math:`0\le k<L`, where :math:`L=\Pi_{n=1}^N l_n`, can be
constructed as

.. math::

  k(i_1,\ldots,i_N) = a_0 + \sum_{n=1}^N a_n i_n

where :math:`a_n`, :math:`i=1,\ldots,N` are coefficients determined by
the particular indexing method. The condition

.. math::

   k(s_1,\ldots,s_N)=0

allows us to determine

.. math::

  a_0 = -\sum_{n=1}^N a_n s_n.

In column-major indexing the region in the :math:`1`-dimensional space
spanned by the first index is contigous:

.. math::

  k(s_1+1,\ldots,i_N) = k(s_1,\ldots,i_N) + 1

Using this condition in the mapping function yields the recursion
relation :math:`a_{j+1}=a_j l_j` with :math:`a_1=1`.

In row-major indexing the region in the :math:`1`-dimensional space
spanned by the last index is contiguous:

.. math::

  k(i_1,\ldots,s_N+1) = k(i_1,\ldots,s_N) + 1

Using this condition in the mapping function yields the recursion
relation :math:`a_{j-1}=a_j l_j` with :math:`a_N=1`.