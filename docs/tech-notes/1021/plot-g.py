from pylab import *
import scipy.special

sz = 24
rc('text', usetex=True)
font = {'size' : sz}
rc('font', **font)

def radau(k, X):
    t1 = scipy.special.eval_legendre(k, X)
    t2 = scipy.special.eval_legendre(k-1, X)
    return 0.5*(-1)**k*(t1-t2)

def gDg(k,X):
    return radau(k,X)

def gGa(k,X):
    return k/(2.*k-1)*radau(k,X) + (k-1)/(2.*k-1)*radau(k-1,X)

def gLo(k,X):
    return (k-1)/(2.*k-1)*radau(k,X) + k/(2.*k-1)*radau(k-1,X)

# make plots k=2
figure(1)
k = 2
X = linspace(-1, 1, 100)
plot(X, gDg(k, X), '-r', label='$g_{DG}$')
plot(X, gGa(k, X), '-k', label='$g_{Ga}$')
plot(X, gLo(k, X), '-m', label='$g_{Lo}$')
tick_params(axis='x', labelsize=sz)
tick_params(axis='y', labelsize=sz)
legend()
xlabel('$\eta$')
ylabel('$g(\eta)$')
title('Correction Polynomial $K=2$')
savefig('gK2.png')

# make plots k=2
figure(2)
k = 4
X = linspace(-1, 1, 100)
plot(X, gDg(k, X), '-r', label='$g_{DG}$')
plot(X, gGa(k, X), '-k', label='$g_{Ga}$')
plot(X, gLo(k, X), '-m', label='$g_{Lo}$')
tick_params(axis='x', labelsize=sz)
tick_params(axis='y', labelsize=sz)
legend()
xlabel('$\eta$')
ylabel('$g(\eta)$')
title('Correction Polynomial $K=4$')
savefig('gK4.png')

show()
