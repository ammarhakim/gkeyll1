import pylab

class Plot2D:
    r"""Plot2D(gkeh : GkeHistoryData, [comp : int, save : bool]) -> Plot2D

    Given a data object, plot it and optionally save it (if ``save``
    is ``True``) to a PNG file.
    """

    def __init__(self, gd, component=0, title=None, save=False, transformMod=None,
                 transformVar=None):

        if transformMod:
            # transform data if needed
            data = transformMod.transformRegistry[transformVar](gd.q)
            defaultTitle = transformVar
        else:
            data = gd.q[:,:,component]
            defaultTitle = gd.base
            
        X = pylab.linspace(gd.lowerBounds[0], gd.upperBounds[0], gd.cells[0]+1)
        Y = pylab.linspace(gd.lowerBounds[1], gd.upperBounds[1], gd.cells[1]+1)
        XX, YY = pylab.meshgrid(X, Y)
        pylab.pcolormesh(XX, YY, data.transpose())
        pylab.colorbar()
        pylab.axis('image')

        if title:
            titleStr = title
        else:
          titleStr = defaultTitle

        extraNm = ""
        if transformMod:
            extraNm = "_" + transformVar
          
        pylab.title('%s' % (titleStr))
        pylab.xlabel('X')
        pylab.ylabel('Y')
        if save:
            figNm = gd.base + extraNm + ("_f%d_c%d.png" % (gd.frame, component))
            pylab.savefig(figNm)

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
          titleStr = gkeh.base    
      pylab.title('History %s[%d]' % (titleStr, component))
      pylab.xlabel('Time')
      pylab.ylabel('%s[%d]' % (titleStr, component))
      if save:
          figNm = gkeh.base + ("_%d.png" % component)
          pylab.savefig(figNm)
