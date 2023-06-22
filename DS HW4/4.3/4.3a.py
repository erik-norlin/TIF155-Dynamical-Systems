# Exercise 4.3a

import numpy as np
import matplotlib.pyplot as plt
import sys

title = '4.3a'

a = 1.4
b = 0.3

cut_tail = 100

T = 1000
dt = 0.01
t = np.arange(0,T,dt)

x = np.zeros_like(t,dtype=float)
y = x.copy()


x[0] = (np.random.uniform()-0.5)
y[0] = (np.random.uniform()-0.5)

for t in range(len(x)-1):
    x[t+1] = y[t]+1-a*x[t]**2
    y[t+1] = b*x[t]

plt.figure(figsize=(7,7))
plt.plot(x[cut_tail:],y[cut_tail:],'.',markersize=0.2)
plt.title(title)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.savefig('Dynamical Systems/DS HW4/4.3/'+title+'.png')
plt.show()