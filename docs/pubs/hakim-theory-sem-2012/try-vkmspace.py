import pylab
import numpy
import math

params = {\
    'text.fontsize': 20,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,    
    }

pylab.rcParams.update(params)

def P0(x):
    return 1

def P1(x):
    return x

def P2(x):
    return 0.5*(3*x*x-1)

# eta -> x
def xeta(e, xc, dx):
    return 0.5*dx*e+xc

def etax(x, xc, dx):
    return 2.0*(x-xc)/dx

def plotLines(X, q0, q1, color):
    dx2 = math.fabs(X[1]-X[0])/2.0
    qlqr = []
    for i in range(X.shape[0]-1):
        ql = q0[i] - q1[i]
        qr = q0[i] + q1[i]
        pylab.plot([X[i], X[i+1]], [ql, qr], color, linewidth=2)
        qlqr.append(ql)
        qlqr.append(qr)
        
    qlqr.sort()
    low, high = qlqr[0], qlqr[len(qlqr)-1]
    for i in range(X.shape[0]-1):
        xe = X[i+1]
        pylab.plot([xe, xe], [low, high], 'm-')

def plotQuads(X, q0, q1, q2, color):
    dx = math.fabs(X[1]-X[0])
    qlqr = []
    for i in range(X.shape[0]-1):
        xc = 0.5*(X[i]+X[i+1])

        ns = 20
        Xs = pylab.linspace(X[i], X[i+1], ns)
        fs = 0.0*Xs
        for m in range(ns):
            et = etax(Xs[m], xc, dx)
            fs[m] = q0[i] + q1[i]*P1(et) + q2[i]*P2(et)

        pylab.plot(Xs, fs, color, linewidth=2)

        qlqr.append(fs.max())
        qlqr.append(fs.min())
        
    qlqr.sort()
    low, high = qlqr[0], qlqr[len(qlqr)-1]
    for i in range(X.shape[0]-1):
        xe = X[i+1]
        pylab.plot([xe, xe], [low, high], 'm-')        


def projectOnLinear(fx, X):
    # weights and absicca for quadrature
    w = [1.0, 1.0]
    eta = [-1/math.sqrt(3.0), 1/math.sqrt(3.0)]
    
    NX = X.shape[0]-1
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    dx = math.fabs(X[1]-X[0])
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = w[0]*P0(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P0(eta[1])*fx(xeta(eta[1], xc, dx))
        f1[i] = w[0]*P1(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P1(eta[1])*fx(xeta(eta[1], xc, dx))

    return 0.5*f0, 1.5*f1

def projectOnQuad(fx, X):
    # weights and absicca for quadrature
    w = [5.0/9.0, 8.0/9.0, 5.0/9.0]
    eta = [-math.sqrt(3.0/5.0), 0, math.sqrt(3.0/5.0)]
    
    NX = X.shape[0]-1
    dx = 0.5*math.fabs(X[1]-X[0])
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    f2 = numpy.zeros((NX,), numpy.float)
    dx = math.fabs(X[1]-X[0])
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = w[0]*P0(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P0(eta[1])*fx(xeta(eta[1], xc, dx)) + w[2]*P0(eta[2])*fx(xeta(eta[2], xc, dx))
        f1[i] = w[0]*P1(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P1(eta[1])*fx(xeta(eta[1], xc, dx)) + w[2]*P1(eta[2])*fx(xeta(eta[2], xc, dx))
        f2[i] = w[0]*P2(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P2(eta[1])*fx(xeta(eta[1], xc, dx)) + w[2]*P2(eta[2])*fx(xeta(eta[2], xc, dx))

    return 0.5*f0, 1.5*f1, 2.5*f2

def initWithLinear(fx, X):
    NX = X.shape[0]-1
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    dx = math.fabs(X[1]-X[0])
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = 0.5*(fx(xeta(1, xc, dx)) + fx(xeta(-1, xc, dx)))
        f1[i] = 0.5*(fx(xeta(1, xc, dx)) - fx(xeta(-1, xc, dx)))

    return f0, f1

def func(x):
    return x**4 + numpy.sin(5*x)

#def func(x):
#    return x**2

X = pylab.linspace(-1, 1, 5+1)
f0p, f1p = projectOnLinear(func, X)
f0q, f1q, f2q = projectOnQuad(func, X)
Xhr = pylab.linspace(-1, 1, 201)
fhr = func(Xhr)

pylab.figure(1)
plotLines(X, f0p, f1p, '-k')
#pylab.plot(Xhr, fhr, '-b')
pylab.title('A Piecewise Linear Function', fontsize='xx-large')
pylab.axis('tight')
pylab.savefig('v1m1.png')

pylab.figure(3)
plotQuads(X, f0q, f1q, f2q, '-k')
#pylab.plot(Xhr, fhr, '-b')
pylab.title('A Piecewise Quadratic Function', fontsize='xx-large')
pylab.axis('tight')
pylab.savefig('v2m1.png')

pylab.show()
