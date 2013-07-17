#!/usr/bin/env python
r"""Command line tool to plot Gkeyll results.

"""

import argparse

# create option parser
parser = argparse.ArgumentParser(description="Plot HDF5 data produced by Gkeyll")
parser.add_argument("files", nargs='+',
                    help="List of HDF5 files to plot. Separate multiple files by spaces.")
parser.add_argument('-c', '--component', default=0, type=int,
                    help='Component of vector data to plot.')
parser.add_argument('-d', '--dg-polyorder', default=0, type=int,
                    help='Polynomial order for DG basis functions.')
parser.add_argument('-g', '--cg-polyorder', default=0, type=int,
                    help='Polynomial order for CG basis functions.')
parser.add_argument('--basis', default='serendip', choices=['serendip', 'serendipity', 
                                                            'tensor-lobatto', 'tensor-gaussian', 
                                                            'tensor-uniform'],
                    help='Type of basis functions used in DG/CG scheme.')
parser.add_argument('--trans-file',
                    help='Name of Python file that performs data transformation.')
parser.add_argument('--trans-var',
                    help='Name of transform to apply from specified transforms file.')
parser.add_argument('--list-trans-vars', action='store_true',
                    help='List all transforms in specified transforms file.')
parser.add_argument('--save-h5', action='store_true',
                    help='Save processed data as HDF5 file.')

args = parser.parse_args()
print args
