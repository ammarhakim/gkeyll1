import pylab
import numpy
import math

# number of cells
NX = 32

def P1(x):
    return x

# weights and absicca for quadrature
w = [1, 1]
eta = [-1/math.sqrt(3), 1/math.sqrt(3)]

# function to differentiate
def fx(x):
    return pylab.sin(x)

def D1fx(x):
    return pylab.cos(x)

def D2fx(x):
    return -pylab.sin(x)

# eta -> x
def xeta(e, xc, dx):
    return 0.5*dx*e+xc

def plotLines(X, q0, q1, color):
    dx = math.fabs(X[1]-X[0])
    for i in range(X.shape[0]-1):
        ql = q0[i] - q1[i]
        qr = q0[i] + q1[i]
        pylab.plot([X[i], X[i+1]], [ql, qr], color)

dx = 2*math.pi/NX
# nodal coordinates
X = pylab.linspace(0, 2*math.pi, NX+1)

# cell averages and slopes of function
f0 = numpy.zeros((NX,), numpy.float)
f1 = numpy.zeros((NX,), numpy.float)

# cell averages and slopes of first derivatives
w0 = numpy.zeros((NX,), numpy.float)
w1 = numpy.zeros((NX,), numpy.float)

# cell averages and slopes of second derivatives
q0 = numpy.zeros((NX,), numpy.float)
q1 = numpy.zeros((NX,), numpy.float)


# project function on piece-wise linear basis functions
for i in range(NX):
    xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

    ## 
    # Use projection on piecewise linear
    ## 
    #f0[i] = 0.5*( w[0]*fx(xeta(eta[0], xc, dx)) + w[1]*fx(xeta(eta[1], xc, dx)) )
    #f1[i] = 1.5*( w[0]*P1(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P1(eta[1])*fx(xeta(eta[1], xc, dx)) )

    ## 
    # Use nodal initialization
    ##     
    f0[i] = 0.5*(fx(xeta(1, xc, dx)) + fx(xeta(-1, xc, dx)))
    f1[i] = 0.5*(fx(xeta(1, xc, dx)) - fx(xeta(-1, xc, dx)))

# compute derivatives in interior
for J in range(NX):
    # funkyness for periodic BCs
    JP = J+1
    JM = J-1
    if (J==0):
        JM = NX-1
    if (J==NX-1):
        JP = 0

    ##
    # 3-point stencil
    ## 
    #w0[J] = -(f1[J]+f0[J]) + (f1[JM]+f0[JM])
    #w1[J] = -3*(f1[J]-f0[J]) - 3*(f1[JM]+f0[JM])

    ##
    # 5-point stencil
    ##     
    w0[J] = 0.5*(f1[JP]-f0[JP]- 2*f1[J] + f1[JM]+f0[JM])
    w1[J] = 0.5*(3*f1[JP]-3*f0[JP] + 6*f0[J] - 3*f1[JM]-3*f0[JM])
    
# compute derivatives in interior
for J in range(NX):
    # funkyness for periodic BCs
    JP = J+1
    JM = J-1
    if (J==0):
        JM = NX-1
    if (J==NX-1):
        JP = 0

    ##
    # 3-point stencil
    ##        
    #q0[J] = (w1[JP]-w0[JP]) - (w1[J]-w0[J])
    #q1[J] = 3*(w1[JP]-w0[JP]) +3*(w1[J]+w0[J])

    ##
    # 5-point stencil
    ##
    q0[J] = 0.5*(w1[JP]-w0[JP] - 2*w1[J] + w1[JP]+w0[JM])
    q1[J] = 0.5*(3*w1[JP]-3*w0[JP] + 6*w0[J] - 3*w1[JM]-3*w0[JM])

# make plots
pylab.figure(1)
pylab.plot(X, fx(X), '-r')
plotLines(X, f0, f1, '-ko')
pylab.axis('tight')

dx = math.fabs(X[1]-X[0])
dx2 = dx*dx

pylab.figure(2)
pylab.plot(X, -D1fx(X), '-r')
plotLines(X, w0/dx, w1/dx, '-k')
pylab.axis('tight')

pylab.figure(3)
pylab.plot(X, D2fx(X), '-r')
plotLines(X, q0/dx2, q1/dx2, '-k')
pylab.axis('tight')

pylab.show()

