import numpy
import pylab
import os
import tables
import gkedata

def plotFile(fileName, argList):
    if not os.path.exists(fileName):
        return False
    fh = tables.openFile(fileName)
    if 'StructGridField' in fh.root:
        # process gridded data
        data = gkedata.GkeData(fileName)
        ndim = data.ndim
        if argList.dg_polyorder:
            pass # doing DG
        elif argList.cg_polyorder:
            pass # doing CG
        else:
            pass # finite-volume data

    elif 'DataStruct' in fh.root:
        # process history data
        pass

    return True


