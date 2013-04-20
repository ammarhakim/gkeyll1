#!/usr/bin/env python
r"""Command line tool to plot Gkeyll results.

"""

import os
import pylab
import gkedata
import gkeplotters

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

# plot history
if options.history:
    hist = gkedata.GkeHistoryData(options.history)
    gkeplotters.PlotHistory(hist, save=options.savePng, title=options.title,
                            component=int(options.component))

# show figure if requested
if options.dontShow == False:
    pylab.show()
