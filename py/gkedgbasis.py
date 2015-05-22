import gkedata
import pylab
import numpy

def makeMatrix(*ll):
    r"""Helper method to contruct a matrix from a list of lists.
    """
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

def interpOnMesh1D(cMat, qIn):
    nInterp, nNodes = cMat.shape[0], cMat.shape[1]
    nx = qIn.shape[0]
    qout = pylab.zeros((nInterp*nx,), pylab.float)
    vList = [qIn[:,i] for i in range(nNodes)]
    for i in range(nInterp):
        qout[i:nInterp*nx:nInterp] = evalSum(cMar[i,:], vList)
    return qout

def intepOnMesh2D(cMat, qIn):
    pass

class GkeDgBasis:
    r"""__init__(dat : GkeData, numNodes : int) -> GkeData

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
        pass


class GkeDgPolyOrder0Basis(GkeDgBasis):
    r"""This is provided to allow treating finite-volume data as DG
    with piecwise contant basis.
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 1)

    def project(self, c):
        return self._getRaw(c)

class GkeDgLobatto1DPolyOrder1Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 1 basis, in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 2)
        self.cMat_i2 = makematrix([0.75,0.25],[0.25,0.75])

    def project(self, c):
        qn = self._getRaw(c)
        return interpOnMesh1D(self.cMat_i2, qn)

class GkeDgLobatto1DPolyOrder2Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 2 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 3)
        self.cMat_i2 = makeMatrix([.5555555555555556,.5555555555555556,-.1111111111111111],[0.0,1.0,0.0],[-.1111111111111111,.5555555555555556,.5555555555555556])

    def project(self, c):
        qn = self._getRaw(c)
        return interpOnMesh1D(self.cMat_i3, qn)

class GkeDgLobatto1DPolyOrder3Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 3 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i4 = makeMatrix([.3964843749999999,.7320061281981992,-.1851311281981993,.05664062500000011],[-0.107421875,.9134865201415714,.2583884798584289,-.06445312499999978],[-.06445312500000017,.2583884798584296,.9134865201415707,-.1074218749999999],[.05664062499999961,-.1851311281981994,.7320061281981994,.3964843749999998])

    def project(self, c):
        qn = self._getRaw(c)
        return interpOnMesh1D(self.cMat_i4, qn)
