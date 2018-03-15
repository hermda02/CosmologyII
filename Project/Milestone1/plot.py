import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True

proj = np.loadtxt('data.dat')
dens = np.loadtxt('densities.dat')

ax = plt.gca()

x_eta = proj[:,0]
a_eta = proj[:,1]
eta   = proj[:,2]
z_eta = proj[:,3]

hx = dens[:,0]
mx = dens[:,1]
bx = dens[:,2]
rx = dens[:,3]
lx = dens[:,4]

#plt.title('Evolution of 'r'$\eta$' )
#plt.plot(a_eta,eta)

#plt.xlabel(r'$a$')
#plt.ylabel(r'$\eta$ ''(cMpc)')

#plt.plot(x_eta,eta)

#plt.xlabel('log'r'$(a)$')
#plt.ylabel(r'$\eta$ ''(cMpc)')

#plt.title('Densities')

#plt.plot(x_eta,mx,label=r"$\Omega_m$")
#plt.plot(x_eta,bx,label=r"$\Omega_b$")
#plt.plot(x_eta,rx,label=r"$\Omega_r$")
#plt.plot(x_eta,lx,label=r"$\Omega_{\Lambda}$")

#plt.xlabel('log'r'$(a)$')
#plt.ylabel('Relative Densities')

#plt.title('Evolution of Hubble Parameter H')
#plt.semilogy(x_eta,hx)

#plt.ylabel('H'r'$(x)$')
#plt.xlabel('log'r'$(a)$')

plt.title('Evolution of Hubble Parameter H')
plt.semilogy(z_eta,hx)

plt.ylabel('H'r'$(z)$')
plt.xlabel(r'$z$')

ax.invert_xaxis()

plt.legend(loc="best")


plt.show()
