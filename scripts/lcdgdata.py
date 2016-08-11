import numpy
import math

class WxRange:
    r"""
    Represents a range of indices.
    """

    def __init__(self, lst):
        self.ndims = len(lst)/2
        self.lower = self.ndims*[0]
        self.upper = self.ndims*[0]
        self.length = self.ndims*[0]        
        for i in range(self.ndims):
            self.lower[i] = lst[2*i]
            self.upper[i] = lst[2*i+1]
            self.length[i] = self.upper[i] - self.lower[i]

class WxColIndexer:
    r"""
    Provides a class to do column-major indexing for arbitrary
    dimensional arrays.
    """

    def __init__(self, r):
        self.ndims = r.ndims
        # allocate memory for coefficients
        self.ai = (self.ndims+1)*[0]

        # set a_1, ..., a_N
        self.ai[1] = 1
        for i in range(2,self.ndims+1):
            self.ai[i] = self.ai[i-1]*r.length[i-2]

        # set a_0
        sm = 0
        for i in range(1,self.ndims+1):
            sm = sm + self.ai[i]*r.lower[i-1]
        self.ai[0] = -sm

    def index1(self, k1):
        return self.ai[0] + k1

    def index2(self, k1, k2):
        return self.ai[0] + k1 + self.ai[2]*k2

    def index3(self, k1, k2, k3):
        return self.ai[0] + k1 + self.ai[2]*k2 + self.ai[3]*k3;

    def index4(self, k1, k2, k3, k4):
        return self.ai[0] + k1 + self.ai[2]*k2 + self.ai[3]*k3 + self.ai[4]*k4

def Pn(n, x):
    p0 = 1.0
    p1 = x

    if n==0: return p0
    if n==1: return p1
    # initiliaze recurrence
    pn = 0.0
    pn1 = p1
    pn2 = p0
    for i in range(2,n+1):
        # use recurrence relation to compute P_n
        pn = (x*(2.*i-1.)*pn1 - (i-1.)*pn2)/(1.*i)
        pn2 = pn1
        pn1 = pn
    return pn

def compQuadPoints(n):
    x = []
    d = 2.0/n
    x.append(-1.0+0.5*d)
    for i in range(1,n):
        x.append(x[i-1]+d)
    return x

def lcDgData(arr, lower, upper, meqn):
    # compute spatial order
    so = arr.shape[1]/meqn
    # construct Legendre polynomials at quadrature points

    xo = compQuadPoints(so)
    legpol = numpy.zeros((so, so), numpy.float)
    for m in range(so):
        for cc in range(so):
            legpol[m,cc] = Pn(m, xo[cc])

    # cell center coordinates
    nx = arr.shape[0]
    dx = (upper-lower)/float(nx)
    Xc = numpy.linspace(lower+0.5*dx, upper-0.5*dx, nx)

    # coordinates at quadrature points
    Xq = numpy.zeros((so*nx,), numpy.float)

    for j in range(so):
        for i in range(nx):
            Xq[so*i+j] = Xc[i] + 0.5*dx*xo[j]

    nxso = nx*so
    # allocate memory for interpolated data
    qres = numpy.zeros((nxso,meqn), numpy.float)
    # interpolate solution to each quadrature point
    rng = WxRange( (0,meqn, 0,so) ) # range object
    idx = WxColIndexer(rng)
    for me in range(meqn):
        for cc in range(so):
            qt = numpy.zeros((nx,), numpy.float)
            for m in range(so):
                qt = qt + legpol[m][cc]*arr[:,idx.index2(me,m)]
            qres[cc:nxso:so,me] = qt

    return Xq, qres
