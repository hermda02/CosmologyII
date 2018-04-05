import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True


elecfrac = np.loadtxt('X_e.dat')
elecred  = np.loadtxt('z_X_e.dat')
logne    = np.loadtxt('logn_e.dat')
taudat   = np.loadtxt('tau.dat')
vis      = np.loadtxt('vis.dat')
vis2     = np.loadtxt('vis2.dat')

ax = plt.gca()

x_rec  = elecfrac[:,0]
X_e    = elecfrac[:,1]

z_rec  = elecred[:,0]
#z_X_e  = elecred[:,1]

#x_rec  = logne[:,0]
logn_e = logne[:,1]

#x_rec   = taudat[:,0]
tau     = taudat[:,1]
tau2    = taudat[:,2]
tau22   = taudat[:,3]

g       = vis[:,1]
g2      = vis2[:,2]
g2 = g2/10.
g22     = vis2[:,3]
g22 = g22/300.
#------------------------------------------------------------------------------------------




plt.title("Electron Fraction")
plt.semilogy(z_rec,X_e)
plt.xlabel('x')
plt.ylabel(r'$X_e$')
plt.xlim(0,1800)
ax.invert_xaxis()
plt.show()

# plt.title("Electron Density")
# plt.semilogy(x_rec,logn_e)
# plt.xlabel('x')
# plt.ylabel('n_e')
# plt.show()

plt.title("Optical Depth")
plt.semilogy(x_rec,tau,label=r'$\tau$')
plt.semilogy(x_rec,tau2,label=r"$\tau`$")
plt.semilogy(x_rec,tau22,label=r"$\tau``$")
plt.xlabel(r'x')
plt.ylabel(r'$\tau$')
plt.xlim(-17,0)
plt.legend(loc="best")
plt.show()


plt.title("Visibility Function")
plt.plot(x_rec,g,label=r'$g$')
plt.plot(x_rec,g2,label=r'$g`$')
plt.plot(x_rec,g22,label=r'$g``$')
plt.xlabel(r'x')
plt.xlim(-7.5,-6)
plt.ylim(-4,5.5)
plt.ylabel(r'g(x)')
plt.legend(loc="best")

#plt.semilogy(z_rec,z_X_e)
#plt.ylabel('X_e')
#plt.xlabel('Redshift')
#plt.xlim(0,1800)
#ax.invert_xaxis()


plt.show()
