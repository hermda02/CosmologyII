import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True
ax = plt.gca()

# ------------------------------------------------------------------------------------------------

int1   = np.loadtxt("thetal_squared1.dat")
int2   = np.loadtxt("thetal_squared2.dat")
int3   = np.loadtxt("thetal_squared3.dat")
int4   = np.loadtxt("thetal_squared4.dat")
int5   = np.loadtxt("thetal_squared5.dat")
int6   = np.loadtxt("thetal_squared6.dat")
tra1   = np.loadtxt("transfer1.dat")
tra2   = np.loadtxt("transfer2.dat")
tra3   = np.loadtxt("transfer3.dat")
tra4   = np.loadtxt("transfer4.dat")
tra5   = np.loadtxt("transfer5.dat")
tra6   = np.loadtxt("transfer6.dat")
C_ls   = np.loadtxt("hi_res_C_ls.dat")
loL    = np.loadtxt("COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",skiprows=3)
hiL    = np.loadtxt("COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt",skiprows=3)

k      = int1[:,0]
theta1 = int1[:,1]
theta2 = int2[:,1]
theta3 = int3[:,1]
theta4 = int4[:,1]
theta5 = int5[:,1]
theta6 = int6[:,1]

loLl   = loL[:,0]
loLCl  = loL[:,1]
yerlolo= loL[:,2]
yerhilo= loL[:,3]
hiLl   = hiL[:,0]
hiLCl  = hiL[:,1]
yerrhi = hiL[:,2]

x      = tra1[:,0]
tran1  = tra1[:,1]
tran2  = tra2[:,1]
tran3  = tra3[:,1]
tran4  = tra4[:,1]
tran5  = tra5[:,1]
tran6  = tra6[:,1]

# l      = C_ls[1:,0]
# C_l    = C_ls[1:,1]

# y = 5775/max(C_l)

# plt.errorbar(loLl,loLCl,yerr=[yerlolo,yerhilo],color='grey')
# plt.errorbar(hiLl,hiLCl,yerr=yerrhi,color='grey')
# plt.plot(l,y*C_l)
# plt.title('Angular Power Spectrum')
# plt.xlabel(r'l')
# plt.ylabel(r'l(l+1)C_l/2\pi')
# plt.xlim(0,1200)
# plt.show()

plt.plot(x,tran1,label=r'l=8')
plt.plot(x,tran2,label=r'l=30')
plt.plot(x,tran3,label=r'l=100')
plt.plot(x,tran4,label=r'l=275')
plt.plot(x,tran5,label=r'l=750')
plt.plot(x,tran6,label=r'l=1200')
plt.legend(loc="best")
plt.xlabel('ck/H_0')
plt.ylabel(r'\Theta_l')
plt.title('Transfer Function (\Theta_l)')
plt.show()

plt.plot(x,theta1,label=r'l=8')
plt.plot(x,theta2,label=r'l=30')
plt.plot(x,theta3,label=r'l=100')
plt.plot(x,theta4,label=r'l=275')
plt.plot(x,theta5,label=r'l=750')
plt.plot(x,theta6,label=r'l=1200')
plt.legend(loc="best")
plt.xlabel('ck/H_0')
plt.ylabel(r'\Theta_l^2/k')
plt.title(r'\Theta_l^2')
plt.show()

