#!/usr/bin/env python
r"""Command line tool to plot Gkeyll results.

"""

import argparse

# create option parser
parser = argparse.ArgumentParser(description="Plot HDF5 data produced by Gkeyll")

dataGroup = parser.add_argument_group('Data specification')
dataGroup.add_argument("files", nargs='+',
                       help="List of HDF5 files to plot. Separate multiple files by spaces. If plotting history"
                       " specify only the first HDF5 file storing the history.")
dataGroup.add_argument('-c', '--component', default=0, type=int,
                       help='Component of vector data to plot. Defaults to 0.')
dataGroup.add_argument('-d', '--dg-polyorder', type=int,
                       help='Polynomial order for DG basis functions.')
dataGroup.add_argument('-g', '--cg-polyorder', type=int,
                       help='Polynomial order for CG basis functions.')
dataGroup.add_argument('--basis', default='serendip', choices=['serendip', 'serendipity', 
                                                               'tensor-lobatto', 'tensor-gaussian', 
                                                               'tensor-uniform'],
                       help='Type of basis functions used in DG/CG scheme. Defaults to serendip.')
dataGroup.add_argument('--trans-file',
                       help='Name of Python file that performs data transformation.')
dataGroup.add_argument('--trans-var',
                       help='Name of transform to apply from specified transforms file.')
dataGroup.add_argument('--list-trans-vars', action='store_true',
                       help='List all transforms in specified transforms file.')
dataGroup.add_argument('--save-h5', action='store_true',
                       help='Save processed data as HDF5 file.')
dataGroup.add_argument('--save-h5-as',
                       help='Optional name of HDF5 file to write.')
dataGroup.add_argument('--save-history',
                       help='Save history data to text file.')

plotGroup = parser.add_argument_group('Plot specification')
plotGroup.add_argument("--title",
                       help="Title to put on plot")
plotGroup.add_argument("--xlabel",
                       help="X-label on plot")
plotGroup.add_argument("--ylabel",
                       help="Y-label on plot")
plotGroup.add_argument("--semilogx", action='store_true', default=False,
                       help="For 1D plots, make a log plot on X-axis.")
plotGroup.add_argument("--semilogy", action='store_true', default=False,
                       help="For 1D plots, make a log plot on Y-axis.")
plotGroup.add_argument("--loglog", action='store_true', default=False,
                       help="For 1D plots, make a log plot on both axis.")
plotGroup.add_argument("--axis-image", action='store_true', default=True,
                       help="For 2D plots preserve aspect ratio. This is the default option.")
plotGroup.add_argument("--axis-tight", action='store_true', default=False,
                       help="For 2D plots do not preserve aspect ratio.")
plotGroup.add_argument("--no-colorbar", action='store_true', default=False,
                       help="For 2D plots a colorbar is not shown.")
plotGroup.add_argument("--dont-show", action='store_true', default=False,
                       help="If set, the plot is not displayed on screen.")
plotGroup.add_argument('-s', "--save", action='store_true', default=False,
                       help="Save plot as a PNG file.")
plotGroup.add_argument("--save-as",
                       help="Optional name of PNG file to save.")

args = parser.parse_args()

