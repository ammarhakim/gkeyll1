import pylab
import numpy
import gkedata
from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 18
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
#rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

def colorbar_adj(obj, mode=1, redraw=False, _fig_=None, _ax_=None, aspect=None):
    '''
    Add a colorbar adjacent to obj, with a matching height
    For use of aspect, see http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.set_aspect ; E.g., to fill the rectangle, try "auto"
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if mode == 1:
        _fig_ = obj.figure; _ax_ = obj.axes
    elif mode == 2: # assume obj is in the current figure/axis instance
        _fig_ = plt.gcf(); _ax_ = plt.gca()
    _divider_ = make_axes_locatable(_ax_)
    _cax_ = _divider_.append_axes("right", size="5%", pad=0.05)
    _cbar_ =  _fig_.colorbar(obj, cax=_cax_)
    if aspect != None:
        _ax_.set_aspect(aspect)
    if redraw:
        _fig_.canvas.draw()
    return _cbar_


class MakeTitle:
    def __init__(self, gd, component, title, transformMod, transformVar, outNm):
        self.title = gd.fName[:-3]+"["+str(component)+"]"
        self.title = self.title + (" at t %g" % gd.time)
        if transformMod:
            self.title = transformVar
            self.title = self.title + (" at t %g" % gd.time)
        if title:
            self.title = title

        self.figName = gd.fName[:-3]
        if transformMod:
            self.figName = self.figName+"_"+transformVar
        else:
            self.figName = self.figName + "-c" + str(component)
        self.figName = self.figName+".png"

        if outNm:
            self.figName = outNm+".png"

class Plot1D:
    r"""Plot1D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)

        if transformMod:
            data = transformMod.transformRegistry[transformVar](gd.q)
        else:
            data = gd.q[:,component]

        dx = (gd.upperBounds[0]-gd.lowerBounds[0])/gd.cells[0]
        X = pylab.linspace(gd.lowerBounds[0]+0.5*dx, gd.upperBounds[0]-0.5*dx, gd.cells[0])
        pylab.plot(X, data)
        pylab.axis('tight')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        if save:
            pylab.savefig(mtitle.figName)

class Plot2D:
    r"""Plot2D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None, options=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)
        if transformMod:
            # transform data if needed
            data = transformMod.transformRegistry[transformVar](gd.q)
        else:
            data = gd.q[:,:,component]

        maskField = 0.0*data # default mask
        if options.maskField:
            # read maskfield
            maskField = gkedata.GkeData(options.maskField).q[:,:,0]

        data = numpy.ma.masked_where(maskField < 0.0, data)
            
        X = pylab.linspace(gd.lowerBounds[0], gd.upperBounds[0], gd.cells[0]+1)
        Y = pylab.linspace(gd.lowerBounds[1], gd.upperBounds[1], gd.cells[1]+1)
        XX, YY = pylab.meshgrid(X, Y)
        im = pylab.pcolormesh(XX, YY, data.transpose())
        if options.cmap:
            plt.set_cmap(options.cmap)
        if options.axisFree:
            pylab.axis('tight')
        else:
            pylab.axis('image')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        colorbar_adj(im)
        plt.tight_layout()
        if save:
            pylab.savefig(mtitle.figName, bbox_inches='tight')

class PlotDg1D:
    r"""PlotDg1D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, dgOrder=1, projOrder=-1, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)

        # number of nodes
        if dgOrder == 1:
            nNodes = 2
        elif dgOrder == 2:
            nNodes = 3
        elif dgOrder == 3:
            raise Exception("1D plotting not implemented for DG polyOrder 3!")
        # number of equations
        numEqns = gd.q.shape[1]/nNodes

        # create data to be plotted
        rawData = numpy.zeros((gd.q.shape[0], nNodes), numpy.float)
        for n in range(nNodes):
            rawData[:,n] = gd.q[:,component+n*numEqns]

        dx = (gd.upperBounds[:]-gd.lowerBounds[:])/gd.cells[:]
        rX = pylab.linspace(gd.lowerBounds[0]+0.5*dx[0],
                            gd.upperBounds[0]-0.5*dx[0], gd.cells[0])
        
        if dgOrder == 1:
            X, data = self.projectOnFinerGrid_f_p1(rX, rawData)
        elif dgOrder == 2:
            X, data = self.projectOnFinerGrid_f_p2(rX, rawData)
        elif dgOrder == 3:
            raise Exception("1D plotting not implemented for DG polyOrder 3!")

        # store
        self.X = X
        self.data = data
        
        pylab.plot(X, data, '-k')
        pylab.axis('tight')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        if save:
            pylab.savefig(mtitle.figName)            

    def evalSum(self, coeff, fields):
        res = 0.0*fields[0]
        for i in range(len(coeff)):
            res = res + coeff[i]*fields[i]
        return res

    def projectOnFinerGrid_f_p1(self, Xc, q):
        dx = Xc[1]-Xc[0]
        nx = Xc.shape[0]

        # mesh coordinates
        xlo = Xc[0]-0.5*dx
        xup = Xc[-1]+0.5*dx
        dx2 = dx/2.0
        Xn = pylab.linspace(xlo+0.5*dx2, xup-0.5*dx2, 2*nx)

        # data
        qn = pylab.zeros((2*Xc.shape[0],), float)
        vList = [q[:,0], q[:,1]]

        # node 1
        c1 = [0.75, 0.25]
        qn[0:2*nx:2] = self.evalSum(c1, vList)
        # node 2
        c2 = [0.25, 0.75]
        qn[1:2*nx:2] = self.evalSum(c2, vList)

        return Xn, qn

    def projectOnFinerGrid_f_p2(self, Xc, q):
        dx = Xc[1]-Xc[0]
        nx = Xc.shape[0]

        # mesh coordinates
        xlo = Xc[0]-0.5*dx
        xup = Xc[-1]+0.5*dx
        dx2 = dx/3.0
        Xn = pylab.linspace(xlo+0.5*dx2, xup-0.5*dx2, 3*nx)

        # data
        qn = pylab.zeros((3*Xc.shape[0],), float)
        vList = [q[:,0], q[:,1], q[:,2]]

        c1 = [5.0/9.0, 5.0/9.0, -1.0/9.0]
        qn[0:3*nx:3] = self.evalSum(c1, vList)
        # node 2
        c2 = [0.0, 1.0, 0.0]
        qn[1:3*nx:3] = self.evalSum(c2, vList)
        # node 3
        c3 = [-1.0/9.0, 5.0/9.0, 5.0/9.0]
        qn[2:3*nx:3] = self.evalSum(c3, vList)

        return Xn, qn    

class PlotDg2D:
    r"""PlotDg2D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, dgOrder=1, projOrder=-1, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None, options=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)

        # number of nodes
        if dgOrder == 1:
            nNodes = 4
        elif dgOrder == 2:
            nNodes = 8
        elif dgOrder == 3:
            raise Exception("1D plotting not implemented for DG polyOrder 3!")
        # number of equations
        numEqns = gd.q.shape[2]/nNodes

        # create data to be plotted
        rawData = numpy.zeros((gd.q.shape[0], gd.q.shape[1], nNodes), numpy.float)
        for n in range(nNodes):
            rawData[:,:,n] = gd.q[:,:,component+n*numEqns]

        dx = (gd.upperBounds[:]-gd.lowerBounds[:])/gd.cells[:]
        rX = pylab.linspace(gd.lowerBounds[0]+0.5*dx[0],
                            gd.upperBounds[0]-0.5*dx[0], gd.cells[0])
        rY = pylab.linspace(gd.lowerBounds[1]+0.5*dx[1], 
                            gd.upperBounds[1]+0.5*dx[1], gd.cells[1])

        if dgOrder == 1:
            if projOrder == 2:
                XX, YY, data = self.projectOnFinerGrid_f24(rX, rY, rawData)
            elif projOrder == 3:
                XX, YY, data = self.projectOnFinerGrid_f29(rX, rY, rawData)
            elif projOrder == 4:
                XX, YY, data = self.projectOnFinerGrid_f216(rX, rY, rawData)
            else:
                raise Exception("ProjOrder of %d not supported!" % projOrder)
        elif dgOrder == 2:
            XX, YY, data = self.projectOnFinerGrid_f39(rX, rY, rawData)

        # store
        self.XX = XX
        self.YY = YY
        self.data = data
                    
        im = pylab.pcolormesh(XX, YY, data.transpose())
        if options.cmap:
            plt.set_cmap(options.cmap)
        if options.axisFree:
            pylab.axis('tight')
        else:
            pylab.axis('image')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        colorbar_adj(im)
        plt.tight_layout()
        if save:
            pylab.savefig(mtitle.figName, bbox_inches='tight')

    def evalSum(self, coeff, fields):
        res = 0.0*fields[0]
        for i in range(len(coeff)):
            res = res + coeff[i]*fields[i]
        return res

    def projectOnFinerGrid_f24(self, Xc, Yc, q):
        dx = Xc[1]-Xc[0]
        dy = Yc[1]-Yc[0]
        nx = Xc.shape[0]
        ny = Yc.shape[0]

        # mesh coordinates
        Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 2*nx+1) # one more
        Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 2*ny+1) # one more
        XXn, YYn = pylab.meshgrid(Xn, Yn)

        # data
        qn = pylab.zeros((2*Xc.shape[0], 2*Yc.shape[0]), float)

        v1 = q[:,:,0]
        v2 = q[:,:,1]
        v3 = q[:,:,2]
        v4 = q[:,:,3]

        vList = [v1,v2,v3,v4]

        # node 1
        c1 = [0.5625,0.1875,0.0625,0.1875]
        qn[0:2*nx:2, 0:2*ny:2] = self.evalSum(c1, vList)

        # node 2
        c2 = [0.1875,0.5625,0.1875,0.0625]
        qn[1:2*nx:2, 0:2*ny:2] = self.evalSum(c2, vList)

        # node 3
        c3 = [0.1875,0.0625,0.1875,0.5625]
        qn[0:2*nx:2, 1:2*ny:2] = self.evalSum(c3, vList)

        # node 4
        c4 = [0.0625,0.1875,0.5625,0.1875]
        qn[1:2*nx:2, 1:2*ny:2] = self.evalSum(c4, vList)
   
        return XXn, YYn, qn                

    def projectOnFinerGrid_f29(self, Xc, Yc, q):
        dx = Xc[1]-Xc[0]
        dy = Yc[1]-Yc[0]
        nx = Xc.shape[0]
        ny = Yc.shape[0]

        # mesh coordinates
        Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 3*nx+1) # one more
        Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 3*ny+1) # one more
        XXn, YYn = pylab.meshgrid(Xn, Yn)

        # data
        qn = pylab.zeros((3*Xc.shape[0], 3*Yc.shape[0]), float)

        v1 = q[:,:,0]
        v2 = q[:,:,1]
        v3 = q[:,:,2]
        v4 = q[:,:,3]

        vList = [v1,v2,v3,v4]


        # node 1
        c1 = [25.0/36.0,5.0/36.0,1.0/36.0,5.0/36.0]
        qn[0:3*nx:3, 0:3*ny:3] = self.evalSum(c1, vList)

        # node 2
        c2 = [5.0/12.0,5.0/12.0,1.0/12.0,1.0/12.0]
        qn[1:3*nx:3, 0:3*ny:3] = self.evalSum(c2, vList)

        # node 3
        c3 = [5.0/36.0,25.0/36.0,5.0/36.0,1.0/36.0]
        qn[2:3*nx:3, 0:3*ny:3] = self.evalSum(c3, vList)

        # node 4
        c4 = [5.0/12.0,1.0/12.0,1.0/12.0,5.0/12.0]
        qn[0:3*nx:3, 1:3*ny:3] = self.evalSum(c4, vList)

        # node 5
        c5 = [1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0]
        qn[1:3*nx:3, 1:3*ny:3] = self.evalSum(c5, vList)

        # node 6
        c6 = [1.0/12.0,5.0/12.0,5.0/12.0,1.0/12.0]
        qn[2:3*nx:3, 1:3*ny:3] = self.evalSum(c6, vList)

        # node 7
        c7 = [5.0/36.0,1.0/36.0,5.0/36.0,25.0/36.0]
        qn[0:3*nx:3, 2:3*ny:3] = self.evalSum(c7, vList)

        # node 8
        c8 = [1.0/12.0,1.0/12.0,5.0/12.0,5.0/12.0]
        qn[1:3*nx:3, 2:3*ny:3] = self.evalSum(c8, vList)

        # node 9
        c9 = [1.0/36.0,5.0/36.0,25.0/36.0,5.0/36.0]
        qn[2:3*nx:3, 2:3*ny:3] = self.evalSum(c9, vList)
   
        return XXn, YYn, qn

    def projectOnFinerGrid_f216(self, Xc, Yc, q):
        dx = Xc[1]-Xc[0]
        dy = Yc[1]-Yc[0]
        nx = Xc.shape[0]
        ny = Yc.shape[0]

        # mesh coordinates
        Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 4*nx+1) # one more
        Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 4*ny+1) # one more
        XXn, YYn = pylab.meshgrid(Xn, Yn)

        # data
        qn = pylab.zeros((4*Xc.shape[0], 4*Yc.shape[0]), float)

        v1 = q[:,:,0]
        v2 = q[:,:,1]
        v3 = q[:,:,2]
        v4 = q[:,:,3]

        vList = [v1,v2,v3,v4]

        ##### node 1
        c1 = [0.765625,0.109375,0.015625,0.109375]
        qn[0:4*nx:4, 0:4*ny:4] = self.evalSum(c1, vList)

        # node 2
        c2 = [0.546875,0.328125,0.046875,0.078125]
        qn[1:4*nx:4, 0:4*ny:4] = self.evalSum(c2, vList)

        # node 3
        c3 = [0.328125,0.546875,0.078125,0.046875]
        qn[2:4*nx:4, 0:4*ny:4] = self.evalSum(c3, vList)

        # node 4
        c4 = [0.109375,0.765625,0.109375,0.015625]
        qn[3:4*nx:4, 0:4*ny:4] = self.evalSum(c4, vList)

        ###### node 5
        c5 = [0.546875,0.078125,0.046875,0.328125]
        qn[0:4*nx:4, 1:4*ny:4] = self.evalSum(c5, vList)

        # node 6
        c6 = [0.390625,0.234375,0.140625,0.234375]
        qn[1:4*nx:4, 1:4*ny:4] = self.evalSum(c6, vList)

        # node 7
        c7 = [0.234375,0.390625,0.234375,0.140625]
        qn[2:4*nx:4, 1:4*ny:4] = self.evalSum(c7, vList)

        # node 8
        c8 = [0.078125,0.546875,0.328125,0.046875]
        qn[3:4*nx:4, 1:4*ny:4] = self.evalSum(c8, vList)

        ###### node 9
        c9 = [0.328125,0.046875,0.078125,0.546875]
        qn[0:4*nx:4, 2:4*ny:4] = self.evalSum(c9, vList)

        # node 10
        c10 = [0.234375,0.140625,0.234375,0.390625]
        qn[1:4*nx:4, 2:4*ny:4] = self.evalSum(c10, vList)

        # node 11
        c11 = [0.140625,0.234375,0.390625,0.234375]
        qn[2:4*nx:4, 2:4*ny:4] = self.evalSum(c11, vList)

        # node 12
        c12 = [0.046875,0.328125,0.546875,0.078125]
        qn[3:4*nx:4, 2:4*ny:4] = self.evalSum(c12, vList)

        ###### node 13
        c13 = [0.109375,0.015625,0.109375,0.765625]
        qn[0:4*nx:4, 3:4*ny:4] = self.evalSum(c13, vList)

        # node 14
        c14 = [0.078125,0.046875,0.328125,0.546875]
        qn[1:4*nx:4, 3:4*ny:4] = self.evalSum(c14, vList)

        # node 15
        c15 = [0.046875,0.078125,0.546875,0.328125]
        qn[2:4*nx:4, 3:4*ny:4] = self.evalSum(c15, vList)

        # node 16
        c16 = [0.015625,0.109375,0.765625,0.109375]
        qn[3:4*nx:4, 3:4*ny:4] = self.evalSum(c16, vList)
   
        return XXn, YYn, qn                

    def projectOnFinerGrid_f39(self, Xc, Yc, q):
        dx = Xc[1]-Xc[0]
        dy = Yc[1]-Yc[0]
        nx = Xc.shape[0]
        ny = Yc.shape[0]

        # mesh coordinates
        Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 3*nx+1) # one more
        Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 3*ny+1) # one more
        XXn, YYn = pylab.meshgrid(Xn, Yn)

        # data
        qn = pylab.zeros((3*Xc.shape[0], 3*Yc.shape[0]), float)

        v1 = q[:,:,0]
        v2 = q[:,:,1]
        v3 = q[:,:,2]
        v4 = q[:,:,3]
        v5 = q[:,:,4]
        v6 = q[:,:,5]
        v7 = q[:,:,6]
        v8 = q[:,:,7]

        vList = [v1,v2,v3,v4,v5,v6,v7,v8]

        # node 1
        c1 = [.2314814814814815,-.1388888888888889,-.06481481481481481,-.1388888888888889,0.462962962962963,.09259259259259259,.09259259259259259,0.462962962962963]
        qn[0:3*nx:3, 0:3*ny:3] = self.evalSum(c1, vList)

        # node 2
        c2 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.8333333333333334,.2777777777777778,.1666666666666667,.2777777777777778]
        qn[1:3*nx:3, 0:3*ny:3] = self.evalSum(c2, vList)

        # node 3
        c3 = [-.1388888888888889,.2314814814814815,-.1388888888888889,-.06481481481481481,0.462962962962963,0.462962962962963,.09259259259259259,.09259259259259259]
        qn[2:3*nx:3, 0:3*ny:3] = self.evalSum(c3, vList)

        # node 4
        c4 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.2777777777777778,.1666666666666667,.2777777777777778,.8333333333333334]
        qn[0:3*nx:3, 1:3*ny:3] = self.evalSum(c4, vList)

        # node 5
        c5 = [-0.25,-0.25,-0.25,-0.25,0.5,0.5,0.5,0.5]
        qn[1:3*nx:3, 1:3*ny:3] = self.evalSum(c5, vList)

        # node 6
        c6 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.2777777777777778,.8333333333333334,.2777777777777778,.1666666666666667]
        qn[2:3*nx:3, 1:3*ny:3] = self.evalSum(c6, vList)

        # node 7
        c7 = [-.1388888888888889,-.06481481481481481,-.1388888888888889,.2314814814814815,.09259259259259259,.09259259259259259,0.462962962962963,0.462962962962963]
        qn[0:3*nx:3, 2:3*ny:3] = self.evalSum(c7, vList)

        # node 8
        c8 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.1666666666666667,.2777777777777778,.8333333333333334,.2777777777777778]
        qn[1:3*nx:3, 2:3*ny:3] = self.evalSum(c8, vList)

        # node 9
        c9 = [-.06481481481481481,-.1388888888888889,.2314814814814815,-.1388888888888889,.09259259259259259,0.462962962962963,0.462962962962963,.09259259259259259]
        qn[2:3*nx:3, 2:3*ny:3] = self.evalSum(c9, vList)
   
        return XXn, YYn, qn            

class PlotHistory:
    r"""PlotHistory(gkeh : GkeHistoryData, [comp : int, save : bool]) -> PlotHistory

    Given a history object, plot it and optionally save it (if
    ``save`` is ``True``) to a PNG file.
    """
    
    def __init__(self, gkeh, component=0, title=None, save=False):
      pylab.plot(gkeh.time, gkeh.history[:,component])
      pylab.axis('tight')
      if title:
          titleStr = title
      else:
          titleStr = gkeh.base + ("[%d]" % component)
      pylab.title('History %s' % (titleStr))
      pylab.xlabel('Time')
      pylab.ylabel('%s' % (titleStr))
      plt.tight_layout()
      if save:
          figNm = gkeh.base + ("_%d.png" % component)
          pylab.savefig(figNm)
