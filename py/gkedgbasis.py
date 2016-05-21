import gkedata
import gkedginterpdat as gid
import pylab
import numpy
import math

## Below are a set of helper functions used in the DG classes

def makeMatrix(*ll):
    nrow = len(ll)
    ncol = len(ll[0])
    mat = numpy.zeros((nrow,ncol), numpy.float)
    for i in range(nrow):
        mat[i,:] = ll[i]
    return mat

def evalSum(coeff, fields):
    res = 0.0*fields[0]
    for i in range(len(coeff)):
        res = res + coeff[i]*fields[i]
    return res

def makeMesh(nInterp, Xc):
    dx = Xc[1]-Xc[0]
    nx = Xc.shape[0]
    xlo = Xc[0]-0.5*dx
    xup = Xc[-1]+0.5*dx
    dx2 = dx/nInterp
    return pylab.linspace(xlo+0.5*dx2, xup-0.5*dx2, nInterp*nx)

def makeMesh2(nInterp, Xc):
    dx = Xc[1]-Xc[0]
    nx = Xc.shape[0]
    xlo = Xc[0]-0.5*dx
    xup = Xc[-1]+0.5*dx
    dx2 = dx/nInterp
    return pylab.linspace(xlo, xup, nInterp*nx+1)

def interpOnMesh1D(cMat, qIn):
    nInterp, nNodes = cMat.shape[0], cMat.shape[1]
    nx = qIn.shape[0]
    qout = pylab.zeros((nInterp*nx,), numpy.float)
    vList = [qIn[:,i] for i in range(nNodes)]
    for i in range(nInterp):
        qout[i:nInterp*nx:nInterp] = evalSum(cMat[i,:], vList)
    return qout

def interpOnMesh2D(cMat, qIn):
    nInterp, nNodes = int(math.sqrt(cMat.shape[0])), cMat.shape[1]
    nx = qIn.shape[0]
    ny = qIn.shape[1]
    qout = pylab.zeros((nInterp*nx,nInterp*ny), numpy.float)
    vList = [qIn[:,:,i] for i in range(nNodes)]
    n = 0
    for j in range(nInterp):
        for i in range(nInterp):
            qout[i:nInterp*nx:nInterp, j:nInterp*ny:nInterp] = evalSum(cMat[n,:], vList)
            n = n+1
    return qout

class GkeDgBasis:
    r"""__init__(dat : GkeData, numNodes : int) -> GkeDgData

    Base class for post-processing DG data. The derived class should
    set the number of nodes in the element.
    """

    def __init__(self, dat, numNodes):
        self.q = dat.q
        self.numNodes = numNodes
        self.numEqns = dat.q.shape[-1]/numNodes
        self.ndim = dat.upperBounds.shape[0]
        self.dx = (dat.upperBounds[:]-dat.lowerBounds[:])/dat.cells[:]
        self.Xc = []
        for d in range(self.ndim):
            self.Xc.append(
                pylab.linspace(dat.lowerBounds[d]+0.5*self.dx[d], dat.upperBounds[d]-0.5*self.dx[d], dat.cells[d])
                )

    def _evalSum(self, coeff, fields):
        r"""evalSum(coeff: [] float, fields: [] array)
        
        Sum arrays in 'fields' list, weighing them with values in 'coeff'
        list.
        """
        return evalSum(coeff, fields)

    def _getRaw(self, component):
        q = self.q
        numEqns = self.numEqns
        shp = [q.shape[i] for i in range(self.ndim)]
        shp.append(self.numNodes)
        rawData = numpy.zeros(shp, numpy.float)
        for n in range(self.numNodes):
            # THERE MUST BE A BETTER WAY TO DO THIS
            if self.ndim == 1:
                rawData[:,n] = q[:,component+n*numEqns]
            elif self.ndim == 2:
                rawData[:,:,n] = q[:,:,component+n*numEqns]
            elif self.ndim == 3:
                rawData[:,:,:,n] = q[:,:,:,component+n*numEqns]
            elif self.ndim == 4:
                rawData[:,:,:,:,n] = q[:,:,:,:,component+n*numEqns]
            elif self.ndim == 5:
                rawData[:,:,:,:,:,n] = q[:,:,:,:,:,component+n*numEqns]
            elif self.ndim == 6:
                rawData[:,:,:,:,:,:,n] = q[:,:,:,:,:,:,component+n*numEqns]

        return rawData        

    def project(self, c):
        r"""project(c : int)
        """
        return 0, 0

## Some notes: The interpolation coefficients are generated in Maxima
## and cut-paste here. Perhaps it may be easier and more compact to
## actually write the Python code to generate these on the
## fly. However, in higher dimensions it may be pretty slow to
## generate these every time.

#################
class GkeDgPolyOrder0Basis(GkeDgBasis):
    r"""This is provided to allow treating finite-volume data as DG
    with piecwise contant basis.
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 1)

    def project(self, c):
        return self.Xc[0], self._getRaw(c)

#################
class GkeDgLobatto1DPolyOrder1Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 1 basis, in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 2)
        self.cMat_i2 = gid.GkeDgLobatto1DPolyOrder1Basis.cMat_i2

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(2, self.Xc[0]), interpOnMesh1D(self.cMat_i2, qn)

#################
class GkeDgLobatto1DPolyOrder2Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 2 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 3)
        self.cMat_i3 = gid.GkeDgLobatto1DPolyOrder2Basis.cMat_i3

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(3, self.Xc[0]), interpOnMesh1D(self.cMat_i3, qn)

#################
class GkeDgLobatto1DPolyOrder3Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 3 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i4 = gid.GkeDgLobatto1DPolyOrder3Basis.cMat_i4

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(4, self.Xc[0]), interpOnMesh1D(self.cMat_i4, qn)

#################
class GkeDgLobatto1DPolyOrder4Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 4 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 5)
        self.cMat_i5 = gid.GkeDgLobatto1DPolyOrder4Basis.cMat_i5

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(5, self.Xc[0]), interpOnMesh1D(self.cMat_i5, qn)

#################
class GkeDgLobatto2DPolyOrder1Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 1 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i2 = gid.GkeDgLobatto2DPolyOrder1Basis.cMat_i2

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(2, self.Xc[0]), makeMesh2(2, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i2, qn)

#################
class GkeDgLobatto2DPolyOrder2Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 2 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 9)
        self.cMat_i3 = gid.GkeDgLobatto2DPolyOrder2Basis.cMat_i3

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(3, self.Xc[0]), makeMesh2(3, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i3, qn)

#################
class GkeDgLobatto2DPolyOrder3Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 3 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 16)
        self.cMat_i4 = gid.GkeDgLobatto2DPolyOrder3Basis.cMat_i4

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(4, self.Xc[0]), makeMesh2(4, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i4, qn)

#################
class GkeDgLobatto2DPolyOrder4Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 4 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 5*5)
        self.cMat_i5 = gid.GkeDgLobatto2DPolyOrder4Basis.cMat_i5

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(5, self.Xc[0]), makeMesh2(5, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i5, qn)

#################
class GkeDgSerendip2DPolyOrder1Basis(GkeDgBasis):
    r"""Serendipity basis (Hakim layout), polyOrder = 1 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i2 = gid.GkeDgSerendip2DPolyOrder1Basis.cMat_i2

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(2, self.Xc[0]), makeMesh2(2, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i2, qn)

#################
class GkeDgSerendip2DPolyOrder2Basis(GkeDgBasis):
    r"""Serendipity basis (Hakim layout), polyOrder = 2 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 8)
        self.cMat_i3 = gid.GkeDgSerendip2DPolyOrder2Basis.cMat_i3

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(3, self.Xc[0]), makeMesh2(3, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i3, qn)

#################
class GkeDgSerendipNorm2DPolyOrder1Basis(GkeDgBasis):
    r"""Serendipity basis (correct, normal layout), polyOrder = 1 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i2 = gid.GkeDgSerendipNorm2DPolyOrder1Basis.cMat_i2

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(2, self.Xc[0]), makeMesh2(2, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i2, qn)

#################
class GkeDgSerendipNorm2DPolyOrder2Basis(GkeDgBasis):
    r"""Serendipity basis (correct, normal layout), polyOrder = 2 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 8)
        self.cMat_i3 = gid.GkeDgSerendipNorm2DPolyOrder2Basis.cMat_i3

    def project(self, c):
        qn = self._getRaw(c)
        X, Y = makeMesh2(3, self.Xc[0]), makeMesh2(3, self.Xc[1])
        XX, YY = pylab.meshgrid(X, Y)
        return XX, YY, interpOnMesh2D(self.cMat_i3, qn)    
