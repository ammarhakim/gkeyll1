import gkedata
import gkedgbasis
from pylab import *


# Lagrange 2D, lobatto polyOrder 1
d = gkedata.GkeData("test2D_lag-2D-p1.h5")
L2d_p1 = gkedgbasis.GkeDgLobatto2DPolyOrder1Basis(d)
Xc, Yc, c0 = L2d_p1.project(0)

figure(1)
pcolormesh(Xc, Yc, transpose(c0))
axis('image')
title('Lagrange 2D, lobatto polyOrder 1')
savefig('L2d_p1.png')

# Lagrange 2D, lobatto polyOrder 2
d = gkedata.GkeData("test2D_lag-2D-p2.h5")
L2d_p2 = gkedgbasis.GkeDgLobatto2DPolyOrder2Basis(d)
Xc, Yc, c0 = L2d_p2.project(0)

figure(2)
pcolormesh(Xc, Yc, transpose(c0))
axis('image')
title('Lagrange 2D, lobatto polyOrder 2')
savefig('L2d_p2.png')

# Lagrange 2D, lobatto polyOrder 3
d = gkedata.GkeData("test2D_lag-2D-p3.h5")
L2d_p3 = gkedgbasis.GkeDgLobatto2DPolyOrder3Basis(d)
Xc, Yc, c0 = L2d_p3.project(0)

figure(3)
pcolormesh(Xc, Yc, transpose(c0))
axis('image')
title('Lagrange 2D, lobatto polyOrder 3')
savefig('L2d_p3.png')

# Lagrange 2D, lobatto polyOrder 4
d = gkedata.GkeData("test2D_lag-2D-p4.h5")
L2d_p4 = gkedgbasis.GkeDgLobatto2DPolyOrder4Basis(d)
Xc, Yc, c0 = L2d_p4.project(0)

figure(4)
pcolormesh(Xc, Yc, transpose(c0))
axis('image')
title('Lagrange 2D, lobatto polyOrder 4')
savefig('L2d_p4.png')

show()
