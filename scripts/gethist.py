#!/usr/bin/env python
#
# def getHist(baseName=None,histName=None,tostdout=0):
# baseName = base name of TxFluids run
# histName = name of history
# tostdout = whether to output to stdout or return as numpy arrays
# 
# 
# Meant to be used as a command line tool (tostdout=1) or embedded 
# within other scripts (tostdout=0)
#
# Example from command line:
# [laptop] pstoltz% ./gethist.py my Bphi_inlet
# 3.26138769318e-10 0.00614445333996
# 6.52277538635e-10 0.012288900236
# 1.03953414002e-09 0.0195847884199
#
# Example embedded:
#In [2]: import gethist
#
#In [3]: [time,data]=gethist.getHist('my','Bphi_inlet')
#
#In [4]: time
#Out[4]: 
#array([  3.26138769e-10,   6.52277539e-10,   1.03953414e-09, ...,
#         9.99552910e-07,   9.99870319e-07,   1.00000000e-06])
#
#In [5]: data
#Out[5]: 
#array([ 0.00614445,  0.0122889 ,  0.01958479, ...,  0.01797907,
#        0.01199911,  0.00955592])
#

import tables
import numpy
import sys

def getHist(baseName=None,histName=None,tostdout=0):
  if not baseName:
    baseName=sys.argv[1]
  if not histName:
    histName=sys.argv[2]
  fname=baseName+'_1.h5'
  hname=histName
  fh=tables.openFile(fname)
  datastr='myb=fh.root.'+hname+'.data[:,0]'
  exec(datastr)
  timestr='myt=fh.root.'+hname+'.timeMesh[:,0]'
  exec(timestr)
  fh.close()
  j=2
  while j:
    try:
      fname=baseName+'_'+str(j)+'.h5'
      fh=tables.openFile(fname)
      datastr='tempb=fh.root.'+hname+'.data[:,0]'
      exec(datastr)
      timestr='tempt=fh.root.'+hname+'.timeMesh[:,0]'
      exec(timestr)
      myb=numpy.append(myb,tempb)
      myt=numpy.append(myt,tempt)
      fh.close()
      j+=1
    except:
      j=0
  if tostdout:
    for j in range(len(myt)):
      print str(myt[j])+' '+str(myb[j])
    return
  else:
    return [myt,myb]

if __name__ == '__main__':
  getHist(tostdout=1)
