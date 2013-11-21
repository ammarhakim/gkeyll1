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

def eps(z,beta_e_val):
    return (m_i/m_e*beta_e_val*z**2-1)*(1+z*plasmaDisp(z)) - kPerpRho**2

def derivEps(z,beta_e_val):
    return (2*m_i/m_e*beta_e_val*z)*(1+z*plasmaDisp(z)) + (m_i/m_e*beta_e_val*z**2-1)*(plasmaDisp(z)+z*derivPlasmaDisp(z))

# Number of points to calculate damping rate at
nPoints = 1000;
beta_e_list = logspace(-6, -1, nPoints)
freqList = zeros(nPoints);
approxFreqList = zeros(nPoints);

m_i = 1.672621777e-27
m_e = 9.10938188e-31
mu0 = 4e-7*pi
eV = 1.602176565e-19

kPerpRho = 0.2
kPar = 0.5
B = 1
Te0 = 250

n0 = B**2*beta_e_list[0]/(2*mu0*Te0*eV)

tol = 1e-4
vA = B/sqrt(mu0*m_i*n0)
# Initial guess for z0 = omega/(k*sqrt(2)*vTe) using approximate expression for wave frequency
z0 = vA/(sqrt(2*Te0*eV/m_e)*sqrt(1+2/beta_e_list[0]*m_e/m_i*kPerpRho**2))

for index, beta_e in enumerate(beta_e_list):
    z0 = optimize.newton(eps,z0,derivEps,(beta_e,),tol,10000)

    nVal = B**2*beta_e/(2*mu0*Te0*eV)
    vA = B/sqrt(mu0*m_i*nVal)
    freqList[index] = fabs(z0.real*sqrt(2*Te0*eV/m_e)/vA);

    # Compute approximate solution using formula
    approxFreqList[index] = 1/sqrt(1 + 2/beta_e_list[index]*m_e/m_i*kPerpRho**2)

# Build list of simulation-derived parameters

nSim = [1.355e16, 2.71e16, 5.42e16, 1.07e17, 4.336e17, 8.673e17, 1.3e18, 8.673e18, 4.336e19, 4.336e20]
tSim = [1.807e-7, 1.834e-7, 1.781e-7, 1.987e-7, 2.496e-7, 3.068e-7, 3.537e-7, 8.397e-7, 1.795e-6, 4.714e-6]
# Really normalized to k*vA
omegaSim = zeros(len(nSim))
betaSim = zeros(len(nSim))

for index, n_val in enumerate(nSim):
  # Divide by 2 to compare with wave freq
  omegaSim[index] = 0.5*(2*pi/tSim[index])/(kPar*B/sqrt(mu0*nSim[index]*m_i))
  betaSim[index] = nSim[index]*2*mu0*Te0*eV/(B**2)

#rc('text', usetex=True)
plt.semilogx(beta_e_list, freqList,'m-',label='Exact')
plt.semilogx(beta_e_list, approxFreqList,'b-',label='Approx')
plt.semilogx(betaSim, omegaSim,'r-o',label='Simulation')

#plt.xlim(betaSim[0],betaSim[-1])
plt.xlabel(r'$\beta_e$')
plt.ylabel(r'$\omega/(k_\parallel v_A)$')
plt.legend(loc='lower right')
plt.autoscale(enable=True,axis='y',tight=True)
plt.savefig('kawDispersionRelation.pdf')
plt.show()
#plt.close()
