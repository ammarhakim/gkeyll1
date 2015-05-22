import gkedata
import pylab
import numpy

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
    return pylab.linspace(xlo+0.5*dx2, xup-0.5*dx2, nInterp*nx+1)

def interpOnMesh1D(cMat, qIn):
    nInterp, nNodes = cMat.shape[0], cMat.shape[1]
    nx = qIn.shape[0]
    qout = pylab.zeros((nInterp*nx,), numpy.float)
    vList = [qIn[:,i] for i in range(nNodes)]
    for i in range(nInterp):
        qout[i:nInterp*nx:nInterp] = evalSum(cMat[i,:], vList)
    return qout

def interpOnMesh2D(cMat, qIn):
    nInterp, nNodes = cMat.shape[0], cMat.shape[1]
    nx = qIn.shape[0]
    ny = qIn.shape[1]
    qout = pylab.zeros((nInterp*nx,nInterp*ny), numpy.float)
    vList = [qIn[:,:,i] for i in range(nNodes)]
    for j in range(nInterp):
        for i in range(nInterp):
            qout[i:nInterp*nx:nInterp, j:nInterp*ny:nInterp] = evalSum(cMat[i,:], vList)
    return qout

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
        return 0, 0

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
        self.cMat_i2 = makeMatrix([0.75,0.25],[0.25,0.75])

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(2, self.Xc[0]), interpOnMesh1D(self.cMat_i2, qn)

#################
class GkeDgLobatto1DPolyOrder2Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 2 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 3)
        self.cMat_i3 = makeMatrix([.5555555555555556,.5555555555555556,-.1111111111111111],[0.0,1.0,0.0],[-.1111111111111111,.5555555555555556,.5555555555555556])

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(3, self.Xc[0]), interpOnMesh1D(self.cMat_i3, qn)

#################
class GkeDgLobatto1DPolyOrder3Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 3 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i4 = makeMatrix([.3964843749999999,.7320061281981992,-.1851311281981993,.05664062500000011],[-0.107421875,.9134865201415714,.2583884798584289,-.06445312499999978],[-.06445312500000017,.2583884798584296,.9134865201415707,-.1074218749999999],[.05664062499999961,-.1851311281981994,.7320061281981994,.3964843749999998])

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(4, self.Xc[0]), interpOnMesh1D(self.cMat_i4, qn)

#################
class GkeDgLobatto1DPolyOrder4Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 4 basis in 1D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 5)
        self.cMat_i5 = makeMatrix([.2663999999999975,.8553363583762964,-0.177600000000006,.08546364162371375,-.02960000000000157],[-.1316000000000021,.7234924181056768,.5263999999999952,-.1746924181056689,.05639999999999855],[2.775557561562891e-17,1.110223024625157e-16,1.0,-5.551115123125783e-16,-1.110223024625157e-16],[.05640000000000095,-.1746924181056717,.5264000000000018,.7234924181056706,-.1315999999999997],[-.02959999999999938,0.0854636416237109,-.1775999999999993,.8553363583762902,.2663999999999995])

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh(5, self.Xc[0]), interpOnMesh1D(self.cMat_i5, qn)

#################
class GkeDgLobatto2DPolyOrder1Basis(GkeDgBasis):
    r"""Lobatto, polyOrder = 1 basis, in 2D
    """

    def __init__(self, dat):
        GkeDgBasis.__init__(self, dat, 4)
        self.cMat_i2 = makeMatrix([0.5625,0.1875,0.1875,0.0625],[0.1875,0.0625,0.5625,0.1875],[0.1875,0.5625,0.0625,0.1875],[0.0625,0.1875,0.1875,0.5625])

    def project(self, c):
        qn = self._getRaw(c)
        return makeMesh2(2, self.Xc[0]), makeMesh2(2, self.Xc[1]), interpOnMesh2D(self.cMat_i2, qn)
