# Numerical solver for a continious dynamical system with streamplot as background

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

xmin = -2
ymin = -2
xmax = 2
ymax = 2

# Streamplot
no_points = 100
x_points = np.linspace(xmin, xmax, no_points)
y_points = np.linspace(ymin, ymax, no_points)
X, Y = np.meshgrid(x_points, y_points)

dxdt_streamplot = Y*(1-X*Y)
dydt_streamplot = X*(X*Y-1)

# Numerical solver
T = 100
dt = 1
t = np.linspace(0,T,int(T/dt)+1)

x0 = 0.01
y0 = 0.01
IC = [x0,y0]

def dynamical_system(IC, t):
    x = IC[0]
    y = IC[1]
    dxdt_solver = x + y - x**2
    dydt_solver = -x + y + 2*x**2
    return [dxdt_solver, dydt_solver]

solver = odeint(dynamical_system, IC, t)
x = solver[:,0]
y = solver[:,1]

# Cutting initial trajectory tail for plotting limit cycles
# cut_tail = 100
# t_new = t[cut_tail:]-cut_tail*dt
# x = solver[cut_tail:,0]
# y = solver[cut_tail:,1]

# Plotting streamplot and solver
fig, ax = plt.subplots(figsize=(7,7))
ax.streamplot(X, Y, dxdt_streamplot, dydt_streamplot, density = 2)
ax.plot(x, y, '-', color='red', linewidth=2)
ax.plot(x0, y0, '.', color='black', markersize=15, label='') # Initial condition
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(xmin,xmax)
ax.set_ylim(xmin,xmax)
ax.set_box_aspect(1) 

title = ''
location = ''

plt.legend(loc="upper left")
plt.savefig(location+title+'.png')
plt.show()