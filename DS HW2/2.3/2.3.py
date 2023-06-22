import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys

xmin = -0.3
ymin = -0.3
xmax = 0.3
ymax = 0.3

# Streamplot
no_points = 1000
x_points = np.linspace(xmin, xmax, no_points)
y_points = np.linspace(ymin, ymax, no_points)
X, Y = np.meshgrid(x_points, y_points)

mu = 0.01
# dxdt_streamplot = mu*X - 5*Y - X**2
# dydt_streamplot = 5*X + mu*Y + 3*Y**3
dxdt_streamplot = mu*X + Y - X**2
dydt_streamplot = -X + mu*Y + 2*X**2

# Numerical integration
T = 100
t = np.linspace(0,T,T*100)
x = np.zeros(T)
y = x.copy()
fp = np.array([[0,0], [0,0]])
x[0] = fp[0,0] + 0.3
y[0] = fp[0,1] + 0.05
T2 = 1000
t2 = np.linspace(0,T2,T2*100)
x2 = np.zeros(T2)
y2 = x2.copy()
x2[0] = fp[0,0] + 0.01
y2[0] = fp[0,1] - 0.01

def dynamical_system(xy, t):
    mu = 0.01
    x = xy[0]
    y = xy[1]
    # dxdt_integration = mu*x - 5*y - x**2
    # dydt_integration = 5*x + mu*y + 3*y**3
    dxdt_integration = mu*x + y - x**2
    dydt_integration = -x + mu*y + 2*x**2
    return [dxdt_integration, dydt_integration]

x0y0 = [x[0],y[0]]
xy = odeint(dynamical_system, x0y0, t)
x = xy[:,0]
y = xy[:,1]
x0y02 = [x2[0],y2[0]]
xy2 = odeint(dynamical_system, x0y02, t2)
x2 = xy2[:,0]
y2 = xy2[:,1]
    
fig, ax = plt.subplots(figsize=(7,7))
ax.streamplot(X, Y, dxdt_streamplot, dydt_streamplot, density = 2)
ax.plot(x, y, '-', color='red', linewidth=3, label='Trajectory going towards the unstable node')
ax.plot(x2, y2, '-', color='orange', linewidth=1, label='Trajectory going away from the unstable node')
ax.plot(fp[0,0],fp[0,1], '.', color='black', markersize=15, label='Unstable node')
# ax.plot(fp[1,0],fp[1,1], '.', color='magenta', markersize=15, label='Unstable spiral')
ax.set_title('2.3: Dynamical system 2, $\mu = 0.01$; indication of suprcritical bifurcation ($a<0$)')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim(xmin,xmax)
ax.set_ylim(xmin,xmax)
ax.set_box_aspect(1) 

plt.legend(loc="upper left")
plt.savefig('23_d2mu0.01.png', bbox_inches='tight')
plt.show()
