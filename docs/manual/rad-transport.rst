Radiation transport
-------------------

Lucee provides algorithms for the solution of the radiation transport
equation (RTE) in various geometries and under different physical
conditions.

Radiation transport in a homogeneous slab
+++++++++++++++++++++++++++++++++++++++++

Radiation transport in homogeneous slabs (plane parallel geometry)
describes several situations of physical interest in ocean and
atmospheric optics. The fundamental problem, also called
Chandrasekhar's basic problem, is described by the equation

.. math::
  :nowrap:

  \notag
  \begin{align}
    \mu\pfrac{L(\tau,\mu)}{\tau} + L(\tau,\mu)
     =
     \frac{\varpi}{2}\sum_{l=0}^L \beta_l P_l(\mu)
     \int_{-1}^1 P_l(\mu')L(\tau,\mu') d\mu' + Q(\tau,\mu)
  \end{align}

for :math:`\tau_0<\tau<\tau_1` and where

* :math:`L(\tau,\mu)` is the radiance in units of 
  :math:`\textrm{Watt m}^{-2}\textrm{sr}^{-1}\textrm{nm}^{-1}`,
* :math:`\tau` is the optical depth,
* :math:`\tau_0` and :math:`\tau_1` are the optical depths of the top and
  bottom surfaces of the slab,
* :math:`\mu` is the cosine of the polar angle measured with the
  positive Z-axis,
* :math:`\varpi` is the albedo of single scattering,
* :math:`Q(\tau,\mu)` are spatially distributed sources.

Further, :math:`\beta_l`, :math:`l=0,\dots,L` are the expansion
coefficients of the phase function :math:`p(\cos\Theta)`, i.e.

.. math::
  :nowrap:

  \notag
  \begin{align}
  p(\cos\Theta) = \sum_{l=0}^L \beta_l P_l(\cos\Theta),
  \end{align}

with the normalization :math:`\beta_0 = 1`.

The integral in the RTE is replaced by a quadrature scheme. Although
any scheme can be used, it is easier to use Gaussian quadrature
separately in the intervals :math:`\mu\in[-1,0]` and
:math:`\mu\in[0,1]`. The sets of weights and abscissa of such a 2Nth
order scheme are denoted by :math:`\{w_i\}` and :math:`\{\mu_i\}`,
:math:`i=1,\ldots,N` respectively.  In the following

.. math::
  :nowrap:

  \notag
  \begin{align}
    \mvec{\Pi}_l &= [P_l(\mu_i),\ldots,P_l(\mu_N)]^T, \\ 
    \mvec{M} &= \textrm{diag}\{\mu_i,\ldots,\mu_N\}, \\
    \mvec{W} &= \textrm{diag}\{w_i,\ldots,w_N\},
  \end{align}

and :math:`\mvec{I}` denotes a :math:`N\times N` unit matrix.