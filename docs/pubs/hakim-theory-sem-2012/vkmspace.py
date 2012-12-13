import pylab
import numpy
import math

params = {\
    'text.usetex': True,
    'text.fontsize': 20,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,    
    }

pylab.rcParams.update(params)

def P1(x):
    return x

# weights and absicca for quadrature
w = [1, 1]
eta = [-1/math.sqrt(3), 1/math.sqrt(3)]

# eta -> x
def xeta(e, xc, dx):
    return 0.5*dx*e+xc

def plotLines(X, q0, q1, color):
    dx2 = math.fabs(X[1]-X[0])/2.0
    qlqr = []
    for i in range(X.shape[0]-1):
        ql = q0[i] - q1[i]
        qr = q0[i] + q1[i]
        pylab.plot([X[i], X[i+1]], [ql, qr], color, linewidth=4)
        qlqr.append(ql)
        qlqr.append(qr)
        
    qlqr.sort()
    low, high = qlqr[0], qlqr[len(qlqr)-1]
    for i in range(X.shape[0]-1):
        xe = X[i+1]
        pylab.plot([xe, xe], [low, high], 'm-')


def projectOnLinear(fx, X):
    NX = X.shape[0]-1
    f0 = numpy.zeros((NX,), numpy.float)
    f1 = numpy.zeros((NX,), numpy.float)
    dx = math.fabs(X[1]-X[0])
    for i in range(NX):
        xc = 0.5*(X[i]+X[i+1]) # cell center coordinate

        f0[i] = 0.5*( w[0]*fx(xeta(eta[0], xc, dx)) + w[1]*fx(xeta(eta[1], xc, dx)) )
        f1[i] = 1.5*( w[0]*P1(eta[0])*fx(xeta(eta[0], xc, dx)) + w[1]*P1(eta[1])*fx(xeta(eta[1], xc, dx)) )

    return f0, f1

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

X = pylab.linspace(-1, 1, 6)
f0p, f1p = projectOnLinear(func, X)
f0c, f1c = initWithLinear(func, X)

dx = math.fabs(X[1]-X[0])
xjm1 = 0.5*(X[1]+X[2])-0.25*dx
xj = 0.5*(X[2]+X[3])
xjp1 = 0.5*(X[3]+X[4])-0.25*dx


pylab.figure(3)
plotLines(X, f0p, f1p, '-ko')
pylab.text(xjm1, 1.75, "$\mathbf{j-1}$", fontsize='xx-large')
pylab.text(xj, 1.75, "$\mathbf{j}$", fontsize='xx-large')
pylab.text(xjp1, 1.75, "$\mathbf{j+1}$", fontsize='xx-large')

xjh_m = xj-0.2*dx
pylab.text(xjh_m, 0.9, "$f_{hj+1/2}^-$", fontsize='xx-large')
xjh_p = xjp1-0.1*dx
pylab.text(xjh_p, 1.2, "$f_{hj+1/2}^+$", fontsize='xx-large')

pylab.axis('tight')
pylab.savefig('v1m1-anno.png')

pylab.show()
