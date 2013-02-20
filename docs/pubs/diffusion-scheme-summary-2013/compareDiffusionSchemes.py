import matplotlib.pyplot as plt
from numpy import *

deltaX = 1
x = linspace(0,2*pi,100) # k_m*deltaX

lambdaDDG    = (1/deltaX)**2*2*(-4+cos(x)+sqrt(6+4*cos(x)-cos(2*x)))
lambdaSDDG   = (1/deltaX)**2*2*(-8-cos(x)+sqrt(42+40*cos(x)-cos(2*x)))
lambdaLDG_AS = (1/deltaX)**2*(-16-2*cos(x)+sqrt(2)*sqrt(93+68*cos(x)+cos(2*x)))
lambdaLDG_S  = (1/deltaX)**2*(-16-2*cos(x)+2*sqrt(42+40*cos(x)-cos(2*x)))
lambdaIR  = (1/deltaX)**2*((sqrt(75*sin(x)**2+64*cos(x)**2+112*cos(x)+49)+cos(x)-16)/2.0)
# Damped eigenvalues:
lambdaDDG1 = (1/deltaX)**2*2*(-4+cos(x)-sqrt(6+4*cos(x)-cos(2*x)))
lambdaSDDG1 = (1/deltaX)**2*2*(-8-cos(x)-sqrt(42+40*cos(x)-cos(2*x)))
lambdaLDG_AS1 = (1/deltaX)**2*(-16-2*cos(x)-sqrt(2)*sqrt(93+68*cos(x)+cos(2*x)))
lambdaLDG_S1  = (1/deltaX)**2*(-16-2*cos(x)-2*sqrt(42+40*cos(x)-cos(2*x)))
lambdaIR1  = -(1/deltaX)**2*((sqrt(75*sin(x)**2+64*cos(x)**2+112*cos(x)+49)-cos(x)+16)/2.0)
# Actual desired behavior
lambdaActual = -(1/deltaX)**2*x**2

plt.figure(1)
plt.plot(x,lambdaDDG,'b',x,lambdaSDDG,'g',x,lambdaLDG_AS,'r',x,lambdaLDG_S,'c')
plt.plot(x, lambdaIR,'k')
plt.plot(x, lambdaActual,'m', linewidth=2)
plt.plot(x, lambdaIR1, '--k')
plt.plot(x,lambdaDDG1,'b--',x,lambdaSDDG1,'g--',x,lambdaLDG_AS1,'r--',x,lambdaLDG_S1,'c--')
plt.xlabel('k$\Delta$x')
plt.ylabel('$\lambda\Delta$x$^2$')
plt.title('Eigenvalue Comparison')
plt.xlim(0,2*pi)
plt.legend(('DDG','SDDG','ALDG','SLDG', 'Recovery', 'Exact'), loc='lower center')

# Put legend outside plot
#ax = plt.subplot(111)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
#ax.legend(('DDG','SDDG','ALDG','SLDG','Exact', 'Recovery'),
#          loc='upper center', bbox_to_anchor=(0.5, -0.1),fancybox=True,shadow=True,ncol=6)

plt.savefig('compareDiffusionSchemes.pdf')
plt.savefig('compareDiffusionSchemes.png')

plt.figure(2)
plt.plot(x,lambdaDDG,'b',x,lambdaSDDG,'g',x,lambdaLDG_AS,'r',x,lambdaLDG_S,'c')
plt.plot(x,lambdaIR,'k')
plt.plot(x,lambdaActual,'m',linewidth=2)
plt.xlabel('k$\Delta$x')
plt.ylabel('$\lambda\Delta$x$^2$')
plt.title('Eigenvalue Comparison')
plt.xlim(0,pi)
plt.ylim(-12,0)
plt.legend(('DDG','SDDG','ALDG','SLDG','Recovery', 'Exact'), loc='center left')

plt.savefig('compareDiffusionSchemes-pi.pdf')

plt.figure(3)
plt.plot(x,lambdaIR,'k', x, lambdaIR1, 'k--')

plt.show()
plt.close()
