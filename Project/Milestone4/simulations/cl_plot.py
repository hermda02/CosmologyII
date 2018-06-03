import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True
ax = plt.gca()

# ------------------------------------------------------------------------------------------------

C_ls   = np.loadtxt("hi_res_C_ls.dat")

C_lm   = np.loadtxt("c_ls_low_m.dat")
C_lm2  = np.loadtxt("c_ls_hi_m.dat")

C_lb   = np.loadtxt("c_ls_low_b.dat")
C_lb2  = np.loadtxt("c_ls_hi_b.dat")

C_lh   = np.loadtxt("c_ls_low_h.dat")
C_lh2  = np.loadtxt("c_ls_hi_h.dat")

C_lr   = np.loadtxt("c_ls_low_r.dat")
C_lr2  = np.loadtxt("c_ls_hi_r.dat")

loL    = np.loadtxt("COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",skiprows=3)
hiL    = np.loadtxt("COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",skiprows=3)

loLl   = loL[:,0]
loLCl  = loL[:,1]
yerlolo= loL[:,2]
yerhilo= loL[:,3]
hiLl   = hiL[:,0]
hiLCl  = hiL[:,1]
yerrhi = hiL[:,2]

l      = C_ls[6:,0]
C_l    = C_ls[6:,1]
Clm    = C_lm[6:,1]
Clm2   = C_lm2[6:,1]
Clb    = C_lb[6:,1]
Clb2   = C_lb2[6:,1]
Clh    = C_lh[6:,1]
Clh2   = C_lh2[6:,1]
Clr    = C_lr[6:,1]
Clr2   = C_lr2[6:,1]

y  = 5775/max(C_l)
y1 = 5775/max(Clm)
y2 = 5775/max(Clm2)
y3 = 5775/max(Clb)
y4 = 5775/max(Clb2)
y5 = 5775/max(Clh)
y6 = 5775/max(Clh2)
y7 = 5775/max(Clr)
y8 = 5775/max(Clr2)

plt.errorbar(loLl,loLCl,yerr=[yerlolo,yerhilo],color='grey',label='Planck')
plt.errorbar(hiLl,hiLCl,yerr=yerrhi,color='grey')
plt.plot(l,y*C_l,label=r'$\Omega_m = 0.224$')
plt.plot(l,y1*Clm,label=r'$\Omega_m = 0.204$')
plt.plot(l,y2*Clm2,label=r'$\Omega_m = 0.244$')
plt.legend(loc="best")
plt.title('Angular Power Spectrum')
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.show()

plt.errorbar(loLl,loLCl,yerr=[yerlolo,yerhilo],color='grey')
plt.errorbar(hiLl,hiLCl,yerr=yerrhi,color='grey')
plt.plot(l,y*C_l,label=r'$\Omega_b = 0.046$')
plt.plot(l,y3*Clb,label=r'$\Omega_b = 0.042$')
plt.plot(l,y4*Clb2,label=r'$\Omega_b = 0.050$')
plt.legend(loc="best")
plt.title('Angular Power Spectrum')
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.show()

plt.errorbar(loLl,loLCl,yerr=[yerlolo,yerhilo],color='grey')
plt.errorbar(hiLl,hiLCl,yerr=yerrhi,color='grey')
plt.plot(l,y*C_l,label=r'$h = 0.7$')
plt.plot(l,y5*Clh,label=r'$h = 0.65$')
plt.plot(l,y6*Clh2,label=r'$h = 0.75$')
plt.legend(loc="best")
plt.title('Angular Power Spectrum')
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.show()

plt.errorbar(loLl,loLCl,yerr=[yerlolo,yerhilo],color='grey')
plt.errorbar(hiLl,hiLCl,yerr=yerrhi,color='grey')
plt.plot(l,y*C_l,label=r'$\Omega_r = 8.3e-5$')
plt.plot(l,y7*Clr,label=r'$\Omega_r = 7.5e-5$')
plt.plot(l,y8*Clr2,label=r'$\Omega_r = 9.0e-5$')
plt.legend(loc="best")
plt.title('Angular Power Spectrum')
plt.xlabel(r'$l$')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlim(0,1200)
plt.show()