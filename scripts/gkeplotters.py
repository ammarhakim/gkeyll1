import pylab

class MakeTitle:
    def __init__(self, gd, component, title, transformMod, transformVar):
        self.title = gd.fName[:-3]+"["+str(component)+"]"
        if transformMod:
            self.title = transformVar
        if title:
            self.title = title

        self.figName = gd.fName[:-3]
        if transformMod:
            self.figName = self.figName+"_"+transformVar
        self.figName = self.figName+".png"

class Plot1D:
    r"""Plot1D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, component=0, title=None, save=False, transformMod=None,
                 transformVar=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar)

        if transformMod:
            data = transformMod.transformRegistry[transformVar](gd.q)
        else:
            data = gd.q[:,component]

        dx = (gd.upperBounds[0]-gd.lowerBounds[0])/gd.cells[0]
        X = pylab.linspace(gd.lowerBounds[0]+0.5*dx, gd.upperBounds[0]-0.5*dx, gd.cells[0])
        pylab.plot(X, data)
        pylab.axis('tight')

        pylab.title('%s at t %g' % (mtitle.title, gd.time))
        pylab.xlabel('X')
        if save:
            pylab.savefig(mtitle.figName)

class Plot2D:
    r"""Plot2D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, component=0, title=None, save=False, transformMod=None,
                 transformVar=None):

        mtitle = MakeTitle(gd, component, title, transformMod, transformVar)
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
        pylab.axis('image')

        pylab.title('%s at t = %g' % (mtitle.title, gd.time))
        pylab.xlabel('X')
        pylab.ylabel('Y')
        if save:
            pylab.savefig(mtitle.figName)

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
