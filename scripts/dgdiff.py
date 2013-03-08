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

def D3fx(x):
    return -pylab.cos(x)

def D5fx(x):
    return pylab.cos(x)

def D4fx(x):
    return pylab.sin(x)

def D5fx(x):
    return pylab.cos(x)

# eta -> x
def xeta(e, xc, dx):
    return 0.5*dx*e+xc

def plotLines(X, q0, q1, color):
    dx2 = math.fabs(X[1]-X[0])/2.0
    for i in range(X.shape[0]-1):
        ql = q0[i] - q1[i]
        qr = q0[i] + q1[i]
        pylab.plot([X[i], X[i+1]], [ql, qr], color)

dx = 2*math.pi/NX
# nodal coordinates
X = pylab.linspace(0, 2*math.pi, NX+1)
# cell center coordinates
Xc = pylab.linspace(0.5*dx, 2*math.pi-0.5*dx, NX)

def projectOnLinearBasis(NX):
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)    
    
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = 0.5*( w[0]*fx(xeta(eta[0], xc, dx)) + w[1]*fx(xeta(eta[1], xc, dx)) )
        f1[i] = 1.5*( w[0]*P1(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P1(eta[1])*fx(xeta(eta[1], xc, dx)) )

    return f0, f1

def evalOnNodes(NX):
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate
        
        f0[i] = 0.5*(fx(xeta(1, xc, dx)) + fx(xeta(-1, xc, dx)))
        f1[i] = 0.5*(fx(xeta(1, xc, dx)) - fx(xeta(-1, xc, dx)))
        
    return f0, f1

def initWithExact(NX):
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    
    for i in range(NX):
        xj = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = (math.cos((2*xj-dx)/2.0)-math.cos((2*xj+dx)/2.0))/dx
        f1[i] = ((2*math.sin((2*xj+dx)/2)-dx*math.cos((2*xj+dx)/2)-2*math.sin((2*xj-dx)/2)-dx*math.cos((2*xj-dx)/2))/dx)*3.0/dx
        
    return f0, f1

# cell averages and slopes of function
#f0, f1 = projectOnLinearBasis(NX)
#f0, f1 = evalOnNodes(NX)
f0, f1 = initWithExact(NX)

def calcGrad(f0, f1):
    NX = f0.shape[0]

    # cell averages and slopes of gradients
    w0 = numpy.zeros((NX,), numpy.float)
    w1 = numpy.zeros((NX,), numpy.float) 

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
        w0[J] = -(f1[J]+f0[J]) + (f1[JM]+f0[JM])
        w1[J] = -3*(f1[J]-f0[J]) - 3*(f1[JM]+f0[JM])

        ##
        # 5-point stencil
        ##     
        #w0[J] = 0.5*(f1[JP]-f0[JP]- 2*f1[J] + f1[JM]+f0[JM])
        #w1[J] = 0.5*(3*f1[JP]-3*f0[JP] + 6*f0[J] - 3*f1[JM]-3*f0[JM])

    return w0, w1

def calcDiv(w0, w1):
    NX = w0.shape[0]

    # cell averages and slopes of gradients
    q0 = numpy.zeros((NX,), numpy.float)
    q1 = numpy.zeros((NX,), numpy.float)

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

    return q0, q1

def calcD2(f0, f1):
    w0, w1 = calcGrad(f0, f1)
    return calcDiv(w0, w1)

def calcD2Direct(f0, f1):
    NX = f0.shape[0]

    # cell averages and slopes D2
    q0 = numpy.zeros((NX,), numpy.float)
    q1 = numpy.zeros((NX,), numpy.float)

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
        q0[J] = -2*f1[JP]+4*f0[JP]-2*f1[J]-8*f0[J]+4*f1[JM]+4*f0[JM]
        q1[J] = -6*f1[JP]+12*f0[JP]-24*f1[J]-6*f0[J]-6*f1[JM]-6*f0[JM]

    return q0, q1

def calcD2DirectSym3(f0, f1):
    NX = f0.shape[0]

    # cell averages and slopes D2
    q0 = numpy.zeros((NX,), numpy.float)
    q1 = numpy.zeros((NX,), numpy.float)

    # compute derivatives in interior
    for J in range(NX):
        # funkyness for periodic BCs
        JP = J+1
        JM = J-1
        if (J==0):
            JM = NX-1
        if (J==NX-1):
            JP = 0

        q0[J] = -6*f1[JP]+8*f0[JP]-16*f0[J]+8*f0[JM]+6*f1[JM]
        q1[J] = -12*f1[JP]+18*f0[JP]-48*f1[J]-18*f0[JM]-12*f1[JM]

    return 0.5*q0, 0.5*q1

def calcD2InterfaceRecovery(f0, f1):
    NX = f0.shape[0]

    # cell averages and slopes D2
    q0 = numpy.zeros((NX,), numpy.float)
    q1 = numpy.zeros((NX,), numpy.float)

    # compute derivatives in interior
    for J in range(NX):
        # funkyness for periodic BCs
        JP = J+1
        JM = J-1
        if (J==0):
            JM = NX-1
        if (J==NX-1):
            JP = 0

        q0[J] = 9*f0[JP]-18*f0[J]+9*f0[JM] - 5*f1[JP]+5*f1[JM]
        q1[J] = 15*f0[JP]-15*f0[JM] - 7*f1[JP]-46*f1[J]-7*f1[JM]

    return 0.25*q0, 0.25*q1

def calcD2DDG(f0, f1):
    NX = f0.shape[0]

    # cell averages and slopes D2
    q0 = numpy.zeros((NX,), numpy.float)
    q1 = numpy.zeros((NX,), numpy.float)

    # compute derivatives in interior
    for J in range(NX):
        # funkyness for periodic BCs
        JP = J+1
        JM = J-1
        if (J==0):
            JM = NX-1
        if (J==NX-1):
            JP = 0

        q0[J] = 2*f0[JM]-4*f0[J]+2*f0[JP] + f1[JM]-f1[JP] 
        q1[J] = -3*f0[JM]+3*f0[JP] - 12*f1[J]

    return q0, q1

Xhr = pylab.linspace(0, 2*math.pi, 1000+1)
# make plots
pylab.figure(1)
pylab.plot(Xhr, fx(Xhr), '-r')
plotLines(X, f0, f1, '-ko')
pylab.axis('tight')
pylab.title('Function')

dx = math.fabs(X[1]-X[0])
dx2 = dx*dx
dx3 = dx*dx*dx
dx4 = dx*dx*dx*dx

w0, w1 = calcGrad(f0, f1)

#pylab.figure(2)
#pylab.plot(X, -D1fx(X), '-r')
#plotLines(X, w0/dx, w1/dx, '-k')
#pylab.plot(Xc, w0/dx, '-m')
#pylab.axis('tight')
#pylab.title('Gradient')

q0, q1 = calcD2(f0, f1)

#pylab.figure(3)
#pylab.plot(X, D2fx(X), '-r')
#plotLines(X, q0/dx2, q1/dx2, '-k')
#pylab.plot(Xc, q0/dx2, '-m')
#pylab.axis('tight')
#pylab.title('Divergence (5-point)')

q0d, q1d = calcD2Direct(f0, f1)
pylab.figure(4)
pylab.plot(Xhr, D2fx(Xhr), '-r')
plotLines(X, q0d/dx2, q1d/dx2, '-k')
pylab.plot(Xc, q0/dx2, '-m')
pylab.plot(Xc, 2/dx*q1d/dx2, 'ob', Xc, D3fx(Xc), 'or')
pylab.plot(Xhr, 6*D2fx(Xhr)/dx, 'r-')
pylab.axis('tight')
pylab.title('Divergence (3-point)')
pylab.savefig('div3-asym.png')

q0s, q1s = calcD2DirectSym3(f0, f1)
pylab.figure(5)
pylab.plot(Xhr, D2fx(Xhr), '-r')
plotLines(X, q0s/dx2, q1s/dx2, '-k')
pylab.plot(Xc, q0s/dx2, '-m')
pylab.plot(Xc, 2/dx*q1s/dx2, 'ob', Xc, D3fx(Xc), 'or')
pylab.plot(Xhr, 168.0/280.0*D3fx(Xhr), 'r-')
pylab.axis('tight')
pylab.title('Divergence (3-point symmetric)')
pylab.savefig('div3-sym.png')

q0i, q1i = calcD2InterfaceRecovery(f0, f1)
pylab.figure(6)
pylab.plot(Xhr, D2fx(Xhr), '-r')
plotLines(X, q0i/dx2, q1i/dx2, '-k')
pylab.plot(Xc, 2/dx*q1i/dx2, 'ob', Xc, D3fx(Xc), 'or')
pylab.plot(Xhr, D3fx(Xhr), 'r-')
pylab.plot(Xc, q0i/dx2, '-m')
pylab.axis('tight')
pylab.title('Divergence (RDG)')
pylab.savefig('div3-rdg.png')

q0dd, q1dd = calcD2DDG(f0, f1)
pylab.figure(7)
pylab.plot(Xhr, D2fx(Xhr), '-r')
plotLines(X, q0dd/dx2, q1dd/dx2, '-k')
pylab.plot(Xc, 2/dx*q1dd/dx2, 'ob', Xc, D3fx(Xc), 'or')
pylab.plot(Xc, q0dd/dx2, '-m')
pylab.axis('tight')
pylab.title('Divergence (DDG)')
pylab.savefig('div3-ddg.png')
    
pylab.show()

