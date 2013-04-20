#!/usr/bin/env python
r"""Command line tool to plot Gkeyll results.

"""

import os
import sys
import pylab
import gkedata
import gkeplotters
import exceptions

from optparse import OptionParser

# set command line options
parser = OptionParser()
parser.add_option('-b', '--base', action = 'store',
                  dest = 'baseName',
                  help = 'Base name of simulation')
parser.add_option('-f', '--frame', action = 'store',
                  dest = 'frame',
                  help = 'Frame to plot.',
                  default = 0)
parser.add_option('-c', '--component', action = 'store',
                  dest = 'component', default=0,
                  help = 'Component to plot')
parser.add_option('-t', '--title', action = 'store',
                  dest = 'title',
                  help = 'Title to put on plots')
parser.add_option('--transforms-file', action = 'store',
                  dest = 'transformsFile',
                  help=  'File for variable transforms')
parser.add_option('-v', '--transform-variable', action = 'store',
                  dest = 'transformVariable',
                  help=  'Name of transform variable to plot')
parser.add_option('-y', '--history', action = 'store',
                  dest = 'history',
                  help = 'Plot specified history')
parser.add_option('-s', '--save', action = 'store_true',
                  dest = 'savePng',
                  help = 'Save png of plot displayed')
parser.add_option('--dont-show', action = 'store_true',
                  dest = 'dontShow',
                  help = 'Do not show plot',
                  default = False)

(options, args) = parser.parse_args()

transformMod = None
# first load in a transforms module, if specified
if options.transformsFile:
    scriptPath = os.path.dirname(sys.argv[0])
    scriptFile = scriptPath + ("/%s.py" % options.transformsFile)
    if os.path.exists(options.transformsFile+".py"):
        sys.path.append(os.path.abspath('.'))
        transformMod = __import__(options.transformsFile)
    elif os.path.exists(scriptFile):
        p = scriptPath + ("/%s" % options.transformsFile)
        sys.path.append(p)
        transformMod = __import__(options.transformsFile)
    else:
        raise exceptions.RuntimeError(
            "Transforms file %s.py does not exist" % options.transformsFile)

# 1D/2D plots
if options.baseName:
    frame = int(options.frame)
    gd = gkedata.GkeData(options.baseName, frame)

    dims = len(gd.q.shape)-1
    # plot, depending on dimension
    if dims == 1:
        pass
    elif dims == 2:
        gkeplotters.Plot2D(gd, save=options.savePng, title=options.title,
                           component=int(options.component), transformMod=transformMod,
                           transformVar=options.transformVariable)
    elif dims == 3:
        raise "Plotting 3D data is not currently supported"

# plot history
if options.history:
    hist = gkedata.GkeHistoryData(options.history)
    gkeplotters.PlotHistory(hist, save=options.savePng, title=options.title,
                            component=int(options.component))

# show figure if requested
if options.dontShow == False:
    pylab.show()
