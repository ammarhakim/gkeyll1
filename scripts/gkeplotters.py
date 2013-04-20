import pylab

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
