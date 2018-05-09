import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

evolution = np.loadtxt("firstdata.dat")

print(evolution[997,0])
print(evolution[998,0])
print(evolution[999,0])

x    = evolution[:999,0]
deltk1 = evolution[:999,1]
deltk10 = evolution[1000:1999,1]
deltk30 = evolution[2000:2999,1]
deltk50 = evolution[3000:3999,1]
deltk80 = evolution[4000:4999,1]
deltk100 = evolution[5000:5999,1]

vk1 = evolution[:999,2]
vk10 = evolution[1000:1999,2]
vk30 = evolution[2000:2999,2]
vk50 = evolution[3000:3999,2]
vk80 = evolution[4000:4999,2]
vk100 = evolution[5000:5999,2]


plt.semilogy(x,deltk1)
plt.semilogy(x,deltk10)
plt.semilogy(x,deltk30)
plt.semilogy(x,deltk50)
plt.semilogy(x,deltk80)
plt.semilogy(x,deltk100)

plt.show()

plt.plot(x,vk1)
plt.plot(x,vk10)
plt.plot(x,vk30)
plt.plot(x,vk50)
plt.plot(x,vk80)
plt.plot(x,vk100)

plt.show()

