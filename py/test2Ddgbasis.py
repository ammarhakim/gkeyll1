import gkedata
import gkedgbasis
from pylab import *

# reference 1D solution
X = linspace(0, 1, 100)
f1d = exp(-75*(X-0.5)**2)

# Lagrange 1D, lobatto polyOrder 1
d = gkedata.GkeData("test1D_lag-1D-p1.h5")
L1d_p1 = gkedgbasis.GkeDgLobatto1DPolyOrder1Basis(d)
Xc, c0 = L1d_p1.project(0)

figure(1)
plot(Xc, c0, 'k-')
plot(X, f1d, 'r-')
title('Lagrange 1D, lobatto polyOrder 1')
savefig('L1d_p1.png')

# Lagrange 1D, lobatto polyOrder 2
d = gkedata.GkeData("test1D_lag-1D-p2.h5")
L1d_p2 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, c0 = L1d_p2.project(0)

figure(2)
plot(Xc, c0, 'k-')
plot(X, f1d, 'r-')
title('Lagrange 1D, lobatto polyOrder 2')
savefig('L1d_p2.png')

# Lagrange 1D, lobatto polyOrder 3
d = gkedata.GkeData("test1D_lag-1D-p3.h5")
L1d_p3 = gkedgbasis.GkeDgLobatto1DPolyOrder3Basis(d)
Xc, c0 = L1d_p3.project(0)

figure(3)
plot(Xc, c0, 'k-')
plot(X, f1d, 'r-')
title('Lagrange 1D, lobatto polyOrder 3')
savefig('L1d_p3.png')

# Lagrange 1D, lobatto polyOrder 4
d = gkedata.GkeData("test1D_lag-1D-p4.h5")
L1d_p4 = gkedgbasis.GkeDgLobatto1DPolyOrder4Basis(d)
Xc, c0 = L1d_p4.project(0)

figure(4)
plot(Xc, c0, 'k-')
plot(X, f1d, 'r-')
title('Lagrange 1D, lobatto polyOrder 4')
savefig('L1d_p4.png')

show()
