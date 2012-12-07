from numpy import *
from scipy import special
from scipy import optimize
import matplotlib.pyplot as plt
import math
from matplotlib import rc

def plasmaDisp(z):
    return 1j*sqrt(pi)*exp(-z**2)*(1+special.erf(1j*z))

def derivPlasmaDisp(z):
    return -2*(1+z*plasmaDisp(z))

def eps(z,T_ratio,T_ratio1):
    return 1 - 1/(2*T_ratio)*derivPlasmaDisp(z)

def derivEps(z,T_ratio,T_ratio1):
    return 1/(T_ratio)*(plasmaDisp(z) + z*derivPlasmaDisp(z))

def isNaN(num):
    return num != num

# Currently not using findRoot in calculation, but it could be useful in future
def findRoot(z0,T_ratio,tol):
    root = 666;
    zCurr = z0;

    while True:
        zNext = zCurr - eps(zCurr,T_ratio,T_ratio)/derivEps(zCurr,T_ratio,T_ratio)
        if isNaN(zCurr.real) or isNaN(zCurr.imag):
            break
        if abs(zNext-zCurr) < tol:
            root = zNext
            break
        zCurr = zNext
    return root

# Number of points to calculate damping rate at
nPoints = 100;
T_ratio_list = linspace(0.1,2.5,nPoints) # Ti/Te
dampingRates = zeros(nPoints);

# Initial guess for z0
z0 = sqrt(1/(2*T_ratio_list[0]))
tol = 1e-6

for index, T_ratio in enumerate(T_ratio_list):
    #z0 = findRoot(z0,T_ratio,tol)
    z0 = optimize.newton(eps,z0,derivEps,(T_ratio,T_ratio),tol,10000)
    if z0 == 666:
        print 'Did not find root for T_ratio = %05f' % T_ratio
        break
    else:
        dampingRates[index] = -z0.imag;

# Import growth rates from input file
dat = loadtxt('simGrowthRatesTest.txt')
simTempRatios = dat[1:,0]
simDampRates  = dat[1:,1:]
simVPoints    = dat[0,1:]

rc('text', usetex=True)
plt.semilogy(T_ratio_list,dampingRates,'m-',label='Exact')

# Plot the different simulation-derived points at the various velocity
# resolutions
for index, vRes in enumerate(simVPoints):
  # Divide simDampRates by two to get field energy damp rates
  yPoints = simDampRates[:,index]*0.5/(sqrt(2)*1.0*0.5)
  plt.semilogy(simTempRatios,yPoints,marker='o',linestyle='None',label=str(vRes))

# plt.xlim(simTempRatios[0],simTempRatios[-1])
plt.xlabel(r'$T = T_i/T_e$')
plt.ylabel('Normalized Damping Rate ($\gamma/\sqrt{2} v_t k$)')
plt.legend(loc='lower right')
#plt.autoscale(enable=True,axis='y',tight=True)
plt.ylim(0.01,1.0)
plt.savefig('ericIonSoundDampingRates.png')
plt.show()
plt.close()
