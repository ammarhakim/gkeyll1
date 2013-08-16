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
        data = gkedata.GkeData(fileName)
    elif 'DataStruct' in fh.root:
        pass

    return True


