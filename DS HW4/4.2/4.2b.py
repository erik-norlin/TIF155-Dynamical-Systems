import numpy as np
import matplotlib.pyplot as plt
import sys

qmin = -20
qmax = 20
q = np.linspace(qmin,qmax,50)
Dq = (np.log(2**q+1)-q*np.log(3))/((1-q)*np.log(3))

plt.figure()
plt.plot(q,Dq)
plt.title('$4.2b$')
plt.xlabel('$q$')
plt.ylabel('$D_q$')
plt.xlim(qmin,qmax)
plt.savefig('Dynamical Systems/DS HW4/4.2/4.2b2.png')
plt.show()