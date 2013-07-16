import pylab

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
            
        X = pylab.linspace(gd.lowerBounds[0], gd.upperBounds[0], gd.cells[0]+1)
        Y = pylab.linspace(gd.lowerBounds[1], gd.upperBounds[1], gd.cells[1]+1)
        XX, YY = pylab.meshgrid(X, Y)
        pylab.pcolormesh(XX, YY, data.transpose())
        pylab.colorbar()
        if options.axisFree:
            pylab.axis('tight')
        else:
            pylab.axis('image')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        if save:
            pylab.savefig(mtitle.figName)

class PlotDg1D:
    r"""PlotDg1D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, dgOrder=1, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)
        rawData = gd.q
        #if transformMod:
        #    # transform data if needed
        #    rawData = transformMod.transformRegistry[transformVar](gd.q)
        #else:
        #    rawData = gd.q[:,:,component]

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

    def __init__(self, gd, dgOrder=1, component=0, title=None, save=False, transformMod=None,
                 transformVar=None, outNm=None, options=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar, outNm)
        rawData = gd.q
        #if transformMod:
        #    # transform data if needed
        #    rawData = transformMod.transformRegistry[transformVar](gd.q)
        #else:
        #    rawData = gd.q[:,:,component]

        dx = (gd.upperBounds[:]-gd.lowerBounds[:])/gd.cells[:]
        rX = pylab.linspace(gd.lowerBounds[0]+0.5*dx[0],
                            gd.upperBounds[0]-0.5*dx[0], gd.cells[0])
        rY = pylab.linspace(gd.lowerBounds[1]+0.5*dx[1], 
                            gd.upperBounds[1]+0.5*dx[1], gd.cells[1])

        if dgOrder == 1:
            XX, YY, data = self.projectOnFinerGrid_f24(rX, rY, rawData)
        elif dgOrder == 2:
            XX, YY, data = self.projectOnFinerGrid_f39(rX, rY, rawData)

        # store
        self.XX = XX
        self.YY = YY
        self.data = data
                    
        pylab.pcolormesh(XX, YY, data.transpose())
        pylab.colorbar()
        if options.axisFree:
            pylab.axis('tight')
        else:
            pylab.axis('image')

        pylab.title(mtitle.title)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        if save:
            pylab.savefig(mtitle.figName)            

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
      if title:
          titleStr = title
      else:
          titleStr = gkeh.base + ("[%d]" % component)
      pylab.title('History %s' % (titleStr))
      pylab.xlabel('Time')
      pylab.ylabel('%s' % (titleStr))
      if save:
          figNm = gkeh.base + ("_%d.png" % component)
          pylab.savefig(figNm)
