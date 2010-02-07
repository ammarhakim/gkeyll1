Radiation transport in a homogeneous slab
-----------------------------------------

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
* :math:`Q(\tau,\mu)` are spatially distributed sources, and
* :math:`P_l(\mu)` is the Legendre polynomial of order :math:`l`.

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
    \mvec{M} &= \textrm{diag}\thinspace \{\mu_i,\ldots,\mu_N\}, \\
    \mvec{W} &= \textrm{diag}\thinspace \{w_i,\ldots,w_N\},
  \end{align}

and :math:`\mvec{I}` denotes a :math:`N\times N` unit matrix.

To find the homogeneous solutions of the RTE, solutions of the form

.. math::
  :nowrap:

  \notag
  \begin{align}
  L_h(\tau,\pm\mu_i) = \phi(\pm\mu_i,\nu) e^{-\tau/\nu}
  \end{align}

are introduced. Here :math:`\phi(\pm\mu_i,\nu)` and :math:`\nu` are
yet to be determined eigenvectors and eigenvalues of the homogeneous
discrete RTE. Denoting

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{\Phi}_{\pm}(\nu) &=
  [\phi(\pm\mu_1,\nu),\ldots,\phi(\pm\mu_N,\nu)]^T, \\
  \mvec{L}_{h\pm}(\tau) &= 
  [L_h(\tau,\pm\mu_1),\ldots,L_h(\tau,\pm\mu_N)]^T
  \end{align}

and using the assumed form in the homogeneous transfer equation it can
be shown that for :math:`\varpi\ne 1` there are exactly 2N eigenvalues
occurring in :math:`\pm` pairs. Denoting this set of eigenvalues by
:math:`\{\pm\nu_j\}, j=1\ldots,N` the eigenvectors are given by

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{\Phi}_+(\nu_j) &= \frac{1}{2}\mvec{M}^{-1}(\mvec{I}+\nu_j \mvec{E})
  \mvec{X}(\lambda_j) \\
  \mvec{\Phi}_-(\nu_j) &= \frac{1}{2}\mvec{M}^{-1}(\mvec{I}-\nu_j \mvec{E})
  \mvec{X}(\lambda_j),
  \end{align}

where :math:`\lambda_j` and :math:`\mvec{X}(\lambda_j)` are the
eigenvalues and eigenvectors of :math:`\mvec{F}\mvec{E}`, and where

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{F} &= 
  \bigg[
  \mvec{I} - \frac{\varpi}{2}
  \sum_{l=0}^L 
  \beta_l \mvec{\Pi}_l\mvec{\Pi}_l^T \mvec{W}(1-(-1)^l)
  \bigg]
  \mvec{M}^{-1}, \\
  \mvec{E} &= 
  \bigg[
  \mvec{I} - \frac{\varpi}{2}
  \sum_{l=0}^L 
  \beta_l \mvec{\Pi}_l\mvec{\Pi}_l^T \mvec{W}(1+(-1)^l)
  \bigg]
  \mvec{M}^{-1},
  \end{align}

and

.. math::
  :nowrap:

  \notag
  \begin{align}
  \pm\nu_j = \pm \lambda_j^{-1/2}.
  \end{align}

The homogeneous solution can be written as a linear combination of the
eigenvectors, i.e.,

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{L}_{h\pm}(\tau) = \sum_{j=1}^N
  \bigg[
  A_j \mvec{\Phi}_{\pm}(\nu_j) e^{-(\tau-\tau_0)/\nu_j}
  +
  B_j \mvec{\Phi}_{\mp}(\nu_j) e^{-(\tau_1-\tau)/\nu_j}
  \bigg],
  \end{align}

where :math:`A_j` and :math:`B_j` are constants to be determined from
the boundary condition.

A general form of the particular solution applicable to any
:math:`Q(\tau,\mu)` can be found using the Green's function
technique. A Green's function for a homogeneous infinite medium is
constructed.  Expressing this in terms of the eigenvectors of the
homogeneous RTE it can be shown that for any source function a
particular solution is given by

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{L}_{p\pm}(\tau) = 
  \sum_{j=1}^N 
  \frac{1}{N(\nu_j)}
  \bigg[
  \script{A}_j(\tau) \mvec{\Phi}_{\pm}(\nu_j)
  +
  \script{B}_j(\tau) \mvec{\Phi}_{\mp}(\nu_j)
  \bigg],
  \end{align}

where

.. math::
  :nowrap:

  \notag
  \begin{align}
  N(\nu_j) = 
  \mvec{\Phi}_+(\nu_j)^T \mvec{W}\mvec{M} \mvec{\Phi}_+(\nu_j)
  -
  \mvec{\Phi}_-(\nu_j)^T \mvec{W}\mvec{M} \mvec{\Phi}_-(\nu_j)
  \end{align}

.. math::
  :nowrap:

  \notag
  \begin{align}
  \script{A}_j(\tau) &=
  \int_{\tau_0}^\tau 
  \bigg[
  \mvec{Q}_+(x)^T \mvec{W} \mvec{\Phi}_+(\nu_j)
  +
  \mvec{Q}_-(x)^T \mvec{W} \mvec{\Phi}_-(\nu_j)
  \bigg]
  e^{-(\tau-x)/\nu_j}
  dx, \\
  \script{B}_j(\tau) &=
  \int_\tau^{\tau_1}
  \bigg[
  \mvec{Q}_+(x)^T \mvec{W} \mvec{\Phi}_-(\nu_j)
  +
  \mvec{Q}_-(x)^T \mvec{W} \mvec{\Phi}_+(\nu_j)
  \bigg]
  e^{-(x-\tau)/\nu_j}
  dx,
  \end{align}

and

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{Q}_\pm(\tau) &= [ Q(\tau,\pm\mu_1),\ldots,Q(\tau,\pm\mu_N)]^T.
  \end{align}

With the particular solution computed the general solution to the RTE
is

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{L}_\pm(\tau) = \mvec{L}_{p\pm}(\tau) +
  \sum_{j=1}^N
  \bigg[
  A_j \gvec{\Phi}_{\pm}(\nu_j) e^{-(\tau-\tau_0)/\nu_j}
  +
  B_j \gvec{\Phi}_{\mp}(\nu_j) e^{-(\tau_1-\tau)/\nu_j}
  \bigg].
  \end{align}

To complete the solution the constants :math:`A_j` and :math:`B_j`
must be determined using the boundary conditions. For example,
consider

.. math::
  :nowrap:

  \notag
  \begin{align}
  L(\tau_0,+\mu) &= 0, \\
  L(\tau_1,-\mu) &= 0,
  \end{align}

for :math:`\mu\in(0,1]`. Using the general solution gives

.. math::
  :nowrap:

  \notag
  \begin{align}
  \sum_{j=1}^N
  \bigg[
  A_j \gvec{\Phi}_+(\nu_j)
  +
  B_j \gvec{\Phi}_-(\nu_j) e^{-(\tau_1-\tau_0)/\nu_j}
  \bigg]
  = - \mvec{L}_{p+}(\tau_0),
  \end{align}

for  :math:`\tau=\tau_0`, and

.. math::
  :nowrap:

  \notag
  \begin{align}
  \sum_{j=1}^N
  \bigg[
  A_j \gvec{\Phi}_-(\nu_j) e^{-(\tau_1-\tau_0)/\nu_j}
  +
  B_j \gvec{\Phi}_+(\nu_j)
  \bigg]
  = - \mvec{L}_{p-}(\tau_1),
  \end{align}

for :math:`\tau=\tau_1`. This system of linear equations can be
inverted to solve for :math:`A_j` and :math:`B_j`.

In general, for other linear boundary conditions one can construct a
system of linear equations written in matix form

.. math::
  :nowrap:

  \notag
  \begin{align}
  \left[
  \begin{matrix}
    \mvec{M}_{11} & \mvec{M}_{12} \\
    \mvec{M}_{21} & \mvec{M}_{22}
  \end{matrix}
  \right]
  \left[
    \begin{matrix}
      \mvec{A} \\
      \mvec{B}
    \end{matrix}
  \right]
  =
  \left[
    \begin{matrix}
      \mvec{P}_1 \\
      \mvec{P}_2
    \end{matrix}
  \right],
  \end{align}

where

.. math::
  :nowrap:

  \notag
  \begin{align}
  \mvec{A} &= [A_1,\ldots,A_N]^T, \\
  \mvec{B} &= [B_1,\ldots,B_N]^T,
  \end{align}

and :math:`\mvec{M}_{ij}`, :math:`\mvec{P}_i`, :math:`i,j=1,2` are
submatrices which depend on boundary conditions at the two surfaces.