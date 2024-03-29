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
parser.add_option('-p', '--plot', action = 'store',
                  dest = 'fileName',
                  help = 'Hdf5 file to plot')
parser.add_option('-m', '--mask', action = 'store',
                  dest = 'maskField',
                  help = 'Field to use as mask (embedded BC simulations only)')
parser.add_option('-c', '--component', action = 'store',
                  dest = 'component', default=0,
                  help = 'Component to plot')
parser.add_option('-t', '--title', action = 'store',
                  dest = 'title',
                  help = 'Title to put on plots')
parser.add_option('--cmap', action = 'store',
                  dest = 'cmap',
                  help = 'Color map to use for 2D plots')
parser.add_option('--transforms-file', action = 'store',
                  dest = 'transformsFile',
                  help=  'File for variable transforms')
parser.add_option('-v', '--transform-variable', action = 'store',
                  dest = 'transformVariable',
                  help=  'Name of transform variable to plot')
parser.add_option('--list-variables', action = 'store_true',
                  dest = 'listTransformVariables',
                  help=  'List all transform variables in file')
parser.add_option('--transform-params', action = 'store',
                  dest = 'transformParams',
                  help = 'List of comma seperated variables needed in computing transforms')
parser.add_option('-y', '--history', action = 'store',
                  dest = 'history',
                  help = 'Plot specified history')
parser.add_option('-s', '--save', action = 'store_true',
                  dest = 'savePng',
                  help = 'Save png of plot displayed')
parser.add_option('-o', '--output', action = 'store',
                  dest = 'outNm',
                  help = 'When saving figures, use this file name')
parser.add_option('--dont-show', action = 'store_true',
                  dest = 'dontShow',
                  help = 'Do not show plot',
                  default = False)
parser.add_option('--project-on', action = 'store',
                  dest = 'projOrder',
                  help = 'Project DG data after splitting cell in these many sub-cells in each direction'
                  'Only makes sense for plotting output from DG schemes')
parser.add_option('-d', '--nodal-order', action = 'store',
                  dest = 'dgOrder',
                  help = 'Polynomial order of DG scheme. '
                  'Only makes sense for plotting output from DG schemes')
parser.add_option('--write-history', action = 'store_true',
                  dest = 'writeHistory',
                  help = 'Write history data',
                  default = False)
parser.add_option('--axis-free', action = 'store_true',
                  dest = 'axisFree',
                  help = "If set, 2D plots won't have equal axis",
                  default = False)

(options, args) = parser.parse_args()

transformMod = None
# first load in a transforms module, if specified
if options.transformsFile:
    varList = []

    if options.transformParams:
        # values needed for transforms to work
        varList = [(x.split('=')[0], float(x.split('=')[1])) for x in options.transformParams.split(",")]

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

    # add global values to module for transforms to work
    for var,val in varList:
        transformMod.__dict__[var] = val

    if options.listTransformVariables:
        print transformMod.transformRegistry.keys()
        exit(0)

# check if to plot data from DG
plotDg = False
dgOrder = 0
if options.dgOrder:
    plotDg = True
    dgOrder = int(options.dgOrder)

projOrder = dgOrder + 1
if options.projOrder:
    projOrder = int(options.projOrder)
        
# 1D/2D plots
if options.fileName:
    gd = gkedata.GkeData(options.fileName)
    dims = len(gd.q.shape)-1

    if plotDg:
        # plot, depending on dimension
        if dims == 1:
            gkeplotters.PlotDg1D(gd, save=options.savePng, title=options.title, 
                                 dgOrder=dgOrder, projOrder=projOrder,
                                 component=int(options.component), transformMod=transformMod,
                                 transformVar=options.transformVariable,
                                 outNm=options.outNm)
        elif dims == 2:
            gkeplotters.PlotDg2D(gd, save=options.savePng, title=options.title,
                                 dgOrder=dgOrder, projOrder=projOrder,
                                 component=int(options.component), transformMod=transformMod,
                                 transformVar=options.transformVariable,
                                 outNm=options.outNm, options=options)
        elif dims == 3:
            raise exceptions.RuntimeError(
                "Plotting 3D DG data is not currently supported")
    else:
        # plot, depending on dimension
        if dims == 1:
            gkeplotters.Plot1D(gd, save=options.savePng, title=options.title,
                               component=int(options.component), transformMod=transformMod,
                               transformVar=options.transformVariable,
                               outNm=options.outNm)
        elif dims == 2:
            gkeplotters.Plot2D(gd, save=options.savePng, title=options.title,
                               component=int(options.component), transformMod=transformMod,
                               transformVar=options.transformVariable,
                               outNm=options.outNm, options=options)
        elif dims == 3:
            raise exceptions.RuntimeError(
                "Plotting 3D data is not currently supported")

# plot history
if options.history:
    hist = gkedata.GkeHistoryData(options.history)
    gkeplotters.PlotHistory(hist, save=options.savePng, title=options.title,
                            component=int(options.component))

    c = int(options.component)
    if options.writeHistory:
        fl = open(hist.base + ".txt", "w")
        for idx in range(hist.time.shape[0]):
            fl.writelines("%.16g %.16g\n" % (hist.time[idx], hist.history[idx,c]))
        fl.close()

# show figure if requested
if options.dontShow == False:
    pylab.show()
