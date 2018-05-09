import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True
ax = plt.gca()

# -----------------------------------------------------------------------------------------

evolution = np.loadtxt("vandb.dat")
phievo    = np.loadtxt('phi_theta.dat')

x        = evolution[:1499,0]

# Dark Matter perturbations
deltk1   = evolution[:1499,1]
deltk10  = evolution[1500:2999,1]
deltk30  = evolution[3000:4499,1]
deltk50  = evolution[4500:5999,1]
deltk80  = evolution[6000:7499,1]
deltk100 = evolution[7500:8999,1]

# Baryon perturbations
deltbk1   = evolution[:1499,2]
deltbk10  = evolution[1500:2999,2]
deltbk30  = evolution[3000:4499,2]
deltbk50  = evolution[4500:5999,2]
deltbk80  = evolution[6000:7499,2]
deltbk100 = evolution[7500:8999,2]

# Dark Matter velocity
vk1   = evolution[:1499,3]
vk10  = evolution[1500:2999,3]
vk30  = evolution[3000:4499,3]
vk50  = evolution[4500:5999,3]
vk80  = evolution[6000:7499,3]
vk100 = evolution[7500:8999,3]

# Baryon velocity
vbk1   = evolution[:1499,4]
vbk10  = evolution[1500:2999,4]
vbk30  = evolution[3000:4499,4]
vbk50  = evolution[4500:5999,4]
vbk80  = evolution[6000:7499,4]
vbk100 = evolution[7500:8999,4]

# Phi
phik1   = phievo[:1499,1]
phik10  = phievo[1500:2999,1]
phik30  = phievo[3000:4499,1]
phik50  = phievo[4500:5999,1]
phik80  = phievo[6000:7499,1]
phik100 = phievo[7500:8999,1]

# Psi
psik1   = phievo[:1499,2]
psik10  = phievo[1500:2999,2]
psik30  = phievo[3000:4499,2]
psik50  = phievo[4500:5999,2]
psik80  = phievo[6000:7499,2]
psik100 = phievo[7500:8999,2]

# Monopole
theta0k1   = phievo[:1499,3]
theta0k10  = phievo[1500:2999,3]
theta0k30  = phievo[3000:4499,3]
theta0k50  = phievo[4500:5999,3]
theta0k80  = phievo[6000:7499,3]
theta0k100 = phievo[7500:8999,3]

# Dipole
theta1k1   = phievo[:1499,4]
theta1k10  = phievo[1500:2999,4]
theta1k30  = phievo[3000:4499,4]
theta1k50  = phievo[4500:5999,4]
theta1k80  = phievo[6000:7499,4]
theta1k100 = phievo[7500:8999,4]

x_tick = [-20,-15,-10,-7.9,-5,0]

plt.semilogy(x,deltk1,label=r"$kc/H_0 = 0.1$")
plt.semilogy(x,deltk10,label=r"$kc/H_0 = 8.36$")
plt.semilogy(x,deltk30,label=r"$kc/H_0 = 85.90$")
plt.semilogy(x,deltk50,label=r"$kc/H_0 = 245.1$")
plt.semilogy(x,deltk80,label=r"$kc/H_0 = 636.8$")
plt.semilogy(x,deltk100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Dark Matter Perturbations")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\delta$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,deltbk1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,deltbk10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,deltbk30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,deltbk50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,deltbk80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,deltbk100,label=r"$kc/H_0 = 1000.0$")
plt.yscale('log')
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Baryon Perturbations")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\delta_b$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,vk1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,vk10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,vk30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,vk50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,vk80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,vk100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Dark Matter Velocity")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$v$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,vbk1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,vbk10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,vbk30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,vbk50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,vbk80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,vbk100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Baryon Velocity")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$v_b$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,phik1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,phik10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,phik30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,phik50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,phik80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,phik100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title(r"$\Phi$ Evolution")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\Phi$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,psik1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,psik10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,psik30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,psik50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,psik80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,psik100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title(r"$\Psi$ Evolution")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\Psi$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,theta0k1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,theta0k10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,theta0k30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,theta0k50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,theta0k80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,theta0k100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Monopole Evolution")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\Theta_0$",size=20)
plt.legend(loc="best")
plt.show()

plt.plot(x,theta1k1,label=r"$kc/H_0 = 0.1$")
plt.plot(x,theta1k10,label=r"$kc/H_0 = 8.36$")
plt.plot(x,theta1k30,label=r"$kc/H_0 = 85.90$")
plt.plot(x,theta1k50,label=r"$kc/H_0 = 245.1$")
plt.plot(x,theta1k80,label=r"$kc/H_0 = 636.8$")
plt.plot(x,theta1k100,label=r"$kc/H_0 = 1000.0$")
plt.text(-7.9, 10000, "Recombination", size=10)
plt.axvspan(-7.397193821,-6.421947418, color='red', alpha=0.5)
ax.set_xticks(x_tick)
ax.set_xticklabels([-20,-15,-10,r"x$_{eq}$",-5,0])
plt.grid()
plt.title("Dipole Evolution")
plt.xlabel(r'$x$',size=20)
plt.ylabel(r"$\Theta_1$",size=20)
plt.legend(loc="best")
plt.show()