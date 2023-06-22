import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys
import os

xmin = -1
ymin = -1
xmax = 1
ymax = 1

# Numerical integration
T = 100
t = np.linspace(0,T,T*100)

x1 = np.zeros(T)
y1 = x1.copy()
x1[0] = 0.01
y1[0] = 0.01

x2 = x1.copy()
y2 = x1.copy()
x2[0] = 0.4
y2[0] = 0.4

x3 = x1.copy()
y3 = x1.copy()

mu = 0.5
omega = 0.5
nu = 0.5
x3[0] = mu
y3[0] = mu

def dynamical_system1(xy, t):
    x = xy[0]
    y = xy[1]
    dxdt_integration = (1/10)*x-y**3-x*y**2-x**2*y-y-x**3
    dydt_integration = x+(1/10)*y+x*y**2+x**3-y**3-x**2*y
    return [dxdt_integration, dydt_integration]

def dynamical_system2(rtheta, t):
    r = rtheta[0]
    drdt_integration = mu*r-r**3
    dthetadt_integration = omega+nu*r**2
    return [drdt_integration, dthetadt_integration]

x0y01 = [x1[0],y1[0]]
xy1 = odeint(dynamical_system1, x0y01, t)
x1 = xy1[:,0]
y1 = xy1[:,1]

x0y02 = [x2[0],y2[0]]
xy2 = odeint(dynamical_system1, x0y02, t)
x2 = xy2[:,0]
y2 = xy2[:,1]

r0theta0 = [mu**0.5,0]
rtheta = odeint(dynamical_system2, r0theta0, t)
r = rtheta[:,0]
theta = rtheta[:,1]
x3 = r*np.cos(theta)
y3 = r*np.sin(theta)
    
fig, ax = plt.subplots(figsize=(7,7))
ax.plot(x1, y1,'-', color='orange', linewidth=2, label='Trajectory going away from f.p. (sys. 2)')
ax.plot(x2, y2,'-', color='red', linewidth=2, label='Trajectory going towards f.p. (sys. 2)')
ax.plot(x3, y3,'-', color='blue', linewidth=2, label='Limit cycle (sys. 1)')
ax.set_title('$3.2b$')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_box_aspect(1) 

plt.legend(loc="upper left")
plt.savefig('Dynamical systems/DS HW3/3.2/3.2b.png', bbox_inches='tight')
plt.show()