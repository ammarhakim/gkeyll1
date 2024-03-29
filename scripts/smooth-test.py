# test smoothing operators
import random
from pylab import *

def P0(x):
    return 1

def P1(x):
    return x

def P2(x):
    return 0.5*(3*x*x-1)

def xeta(e, xc, dx):
    return 0.5*dx*e+xc

def etax(x, xc, dx):
    return 2.0*(x-xc)/dx

def projectSinOnQuadBasis(X):
    dx = X[1]-X[0]
    nx = X.shape[0]
    fx = zeros((nx,3), float)
    cos, sin = math.cos, math.sin
    
    for i in range(nx):
        xc = X[i]
        fx[i,0] = cos((2*xc-dx)/2.0)-cos((2*xc+dx)/2.0)
        fx[i,1] = (2*sin((2*xc+dx)/2.0)-dx*cos((2*xc+dx)/2.0)-2*sin((2*xc-dx)/2.0)-dx*cos((2*xc-dx)/2.0))/dx
        fx[i,2] = (6*dx*sin((2*xc+dx)/2.0)+(12-dx**2)*cos((2*xc+dx)/2.0)+6*dx*sin((2*xc-dx)/2.0)+(dx**2-12)*cos((2*xc-dx)/2.0))/dx**2

    # rescale
    fx[:,0] = 1.0/dx*fx[:,0]
    fx[:,1] = 3.0/dx*fx[:,1]
    fx[:,2] = 5.0/dx*fx[:,2]

    return fx

def projectLinOnQuadBasis(X):
    dx = X[1]-X[0]
    nx = X.shape[0]
    fx = zeros((nx,3), float)
    
    for i in range(nx):
        xc = X[i]
        fx[i,0] = xc
        fx[i,1] = dx/2.0
        fx[i,2] = 0.0

    return fx

def projectQuadOnQuadBasis(X):
    dx = X[1]-X[0]
    nx = X.shape[0]
    fx = zeros((nx,3), float)
    
    for i in range(nx):
        xc = X[i]
        fx[i,2] = dx**2/6
        fx[i,0] = xc**2 + 0.5*fx[i,2]
        fx[i,1] = dx*xc

    return fx

def plotQuads(X, q0, q1, q2, color):
    dx = math.fabs(X[1]-X[0])
    qlqr = []
    for i in range(X.shape[0]-1):
        xc = 0.5*(X[i]+X[i+1])

        ns = 20
        Xs = linspace(X[i], X[i+1], ns)
        fs = 0.0*Xs
        for m in range(ns):
            et = etax(Xs[m], xc, dx)
            fs[m] = q0[i] + q1[i]*P1(et) + q2[i]*P2(et)

        plot(Xs, fs, color, linewidth=2)

        qlqr.append(fs.max())
        qlqr.append(fs.min())
        
    qlqr.sort()
    low, high = qlqr[0], qlqr[len(qlqr)-1]
    for i in range(X.shape[0]-1):
        xe = X[i+1]
        plot([xe, xe], [low, high], '-', color='#808080')

def smoothQuad(fx, c0=2, beta=1.0):
    nx = fx.shape[0]-2
    gx = 0.0*fx
    alpha = 1.0/5.0*(3*beta-1)

    a = (5*alpha-1)*c0+10*alpha+2
    b = (10*alpha+2)*c0+20*alpha-4
    for i in range(1,nx+1):
        gx[i,0] = (fx[i+1,0]+c0*fx[i,0]+fx[i-1,0])/(c0+2) \
            - (fx[i+1,1]-fx[i-1,1])/12.0 \
            + (c0-2)*(fx[i+1,2]-2*fx[i,2]+fx[i-1,2])/(20*c0+40)

        gx[i,1] = (fx[i+1,0]-fx[i-1,0])/4 \
            - beta*(fx[i+1,1]-2*fx[i,1]+fx[i-1,1])/4.0 \
            + alpha*(fx[i+1,2]-fx[i-1,2])/4

        gx[i,2] = (c0-2)/(4*c0+8)*(fx[i+1,0]-2*fx[i,0]+fx[i-1,0]) \
            - (3*beta-1)*(fx[i+1,1]-fx[i-1,1])/12.0 \
            + (a*fx[i+1,2] + b*fx[i,2] + a*fx[i-1,2])/(20*c0+40)

    return gx

def smoothLin(fx):
    nx = fx.shape[0]-2
    gx = 0.0*fx
    for i in range(1,nx+1):
        gx[i,0] = (fx[i+1,0]+2*fx[i,0]+fx[i-1,0])/4 - (fx[i+1,1]-fx[i-1,1])/12.0
        gx[i,1] = (fx[i+1,0]-fx[i-1,0])/4 - (fx[i+1,1]-2*fx[i,1]+fx[i-1,1])/12.0
    return gx

# parameters for plot
nx = 8
Lx = 2*pi
dx = Lx/nx
X = linspace(-0.5*dx, Lx+0.5*dx, nx+2)
Xedge = linspace(-dx, Lx+dx, nx+3)

### f(x) = sin(x)
figure(1)
fx = projectSinOnQuadBasis(X)
plotQuads(Xedge, fx[:,0], fx[:,1], fx[:,2], 'r')

c0 = 10.0
beta = 3.0/4.0
alpha = 1/5.0*(3*beta-1)
plotQuads(Xedge, fx[:,0], beta*fx[:,1], alpha*fx[:,2], '--g')

Xhr = linspace(Xedge[0], Xedge[-1], 100)
plot(Xhr, sin(Xhr), 'b-')

# smooth it
gx = smoothQuad(fx, c0, beta)
plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
title('Red: DG. Black: Smooth')

axis('tight')

#print ("Area under curve f(x)=sin(x)", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

### f(x) = x**2
fx = projectQuadOnQuadBasis(X)
figure(3)
plotQuads(Xedge, fx[:,0], fx[:,1], fx[:,2], 'r')
plotQuads(Xedge, fx[:,0], beta*fx[:,1], alpha*fx[:,2], 'g--')
Xhr = linspace(Xedge[0], Xedge[-1], 100)
plot(Xhr, Xhr**2, 'b-')

# smooth it
gx = smoothQuad(fx, c0, beta)
plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
title('Red: DG. Black: Smooth')

axis('tight')

#print ("Area under curve f(x)=x^2", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

### f(x) = x
fx = projectLinOnQuadBasis(X)
figure(2)
plotQuads(Xedge, fx[:,0], fx[:,1], fx[:,2], 'r')
#plotQuads(Xedge, fx[:,0], fx[:,1], 2.0/5.0*fx[:,2], 'g--')

# smooth it
gx = smoothQuad(fx)
plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
title('Red: DG. Black: Smooth')

axis('tight')

#print ("Area under curve f(x)=x", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

### f(x) = sin(x) with perturbations
figure(4)
fx = projectSinOnQuadBasis(X)
for i in range(fx.shape[0]):
    fx[i,2] = 2.0*(random()-0.5)*fx[i,2]
plotQuads(Xedge, fx[:,0], fx[:,1], fx[:,2], 'r')

c0 = 10.0
beta = 3.0/4.0
alpha = 1/5.0*(3*beta-1)
plotQuads(Xedge, fx[:,0], beta*fx[:,1], alpha*fx[:,2], '--g')

Xhr = linspace(Xedge[0], Xedge[-1], 100)
plot(Xhr, sin(Xhr), 'b-')

# smooth it
gx = smoothQuad(fx, c0, beta)
plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
title('Red: DG. Black: Smooth')

axis('tight')

#print ("Area under curve f(x)=sin(x)", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

# ######################## Linear tests
# nx = 3
# Lx = 2*pi
# dx = Lx/nx
# X = linspace(-0.5*dx, Lx+0.5*dx, nx+2)
# Xedge = linspace(-dx, Lx+dx, nx+3)

# ### f(x) = sin(x) 
# fx = projectSinOnQuadBasis(X)
# figure(4)
# plotQuads(Xedge, fx[:,0], fx[:,1], 0.0*fx[:,2], 'r')
# plotQuads(Xedge, fx[:,0], 1.0/3.0*fx[:,1], 0.0/5.0*fx[:,2], 'g--')

# Xhr = linspace(Xedge[0], Xedge[-1], 100)
# #plot(Xhr, sin(Xhr), 'b-')

# # smooth it
# gx = smoothLin(fx)
# plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
# title('Linear. Red: DG. Black: Smooth')

# axis('tight')

# print ("Area under curve f(x)=sin(x)", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

# ### f(x) = x**2
# fx = projectQuadOnQuadBasis(X)
# figure(5)
# plotQuads(Xedge, fx[:,0], fx[:,1], 0.0*fx[:,2], 'r')
# plotQuads(Xedge, fx[:,0], 1.0/3.0*fx[:,1], 0.0*fx[:,2], 'g--')

# # smooth it
# gx = smoothLin(fx)
# plotQuads(Xedge[1:-1], gx[1:-1,0], gx[1:-1,1], gx[1:-1,2], 'k')
# title('Linear. Red: DG. Black: Smooth')

# axis('tight')

# print ("Area under curve f(x)=x^2", sum(fx[1:-1,0]), sum(gx[1:-1,0]))

show()
