import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys
import os

xmin = -1
ymin = -1
zmin = -1
xmax = 1
ymax = 1
zmax = 1

# Streamplot
# no_points = 100
# x_points = np.linspace(xmin, xmax, no_points)
# y_points = np.linspace(ymin, ymax, no_points)
# z_points = np.linspace(zmin, zmax, no_points)
# X, Y, Z = np.meshgrid(x_points, y_points)

# dxdt_streamplot = mu*X + Y - X**2
# dydt_streamplot = -X + mu*Y + 2*X**2
# dzdt_streamplot = 

# Numerical integration
T = 500
t = np.linspace(0,T,T*100)
x = np.zeros(T)
y = x.copy()
z = x.copy()
# fp = np.array([[(1+(mu**2))/(2+mu), (-(2*mu-1)*(1+mu**2))/((2+mu)**2)], [0,0]])
# x[0] = fp[0,0] - 0.01
# y[0] = fp[0,1] - 0.01
x[0] = 0.01
y[0] = 0.01
z[0] = 0.01

def dynamical_system(xyz, t):
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    dxdt_integration = 10*(y-x)
    dydt_integration = 28*x-y-x*z
    dzdt_integration =  x*y-(8/3)*z
    return [dxdt_integration, dydt_integration, dzdt_integration]

x0y0z0 = [x[0],y[0],z[0]]
xyz = odeint(dynamical_system, x0y0z0, t)
x = xyz[:,0]
y = xyz[:,1]
z = xyz[:,2]
    
fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(7,7))
# ax.streamplot(X, Y, dxdt_streamplot, dydt_streamplot, density = 2)
ax.plot(x[100:], y[100:], z[100:], '-', color='red', linewidth=0.25)
# ax.plot(fp[0,0],fp[0,1], '.', color='black', markersize=15, label='Saddle node')
# ax.plot(fp[1,0],fp[1,1], '.', color='magenta', markersize=15, label='Unstable spiral')
ax.set_title('$3.1b$')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# ax.set_xlim(xmin,xmax)
# ax.set_ylim(ymin,ymax)
# ax.set_ylim(zmin,zmax)
# ax.set_box_aspect(1) 

# plt.legend(loc="upper left")
script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, '3.1/')
plt.savefig('Dynamical systems/DS HW3/3.1/3.1b.png', bbox_inches='tight')
plt.show()