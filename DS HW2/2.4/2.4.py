import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys


xmin = -0.5
ymin = -0.5
xmax = 0.75
ymax = 0.75

# Streamplot
no_points = 100
x_points = np.linspace(xmin, xmax, no_points)
y_points = np.linspace(ymin, ymax, no_points)
X, Y = np.meshgrid(x_points, y_points)

global mu 
mu= 0.06599
dxdt_streamplot = mu*X + Y - X**2
dydt_streamplot = -X + mu*Y + 2*X**2

# Numerical integration
T = 100
t = np.linspace(0,T,T*10)
x = np.zeros(T)
y = x.copy()
x2 = x.copy()
y2 = y.copy()
fp = np.array([[(1+(mu**2))/(2+mu), (-(2*mu-1)*(1+mu**2))/((2+mu)**2)], [0,0]])
# x[0] = fp[0,0] - 0.01
# y[0] = fp[0,1] - 0.01
x[0] = 0.01
y[0] = 0.01
# x2[0] = fp[1,0] - 0.005
# y2[0] = fp[1,1] - 0.005

def dynamical_system(xy, t):
    # mu = 0.065
    x = xy[0]
    y = xy[1]
    dxdt_integration = mu*x + y - x**2
    dydt_integration = -x + mu*y + 2*x**2
    return [dxdt_integration, dydt_integration]

x0y0 = [x[0],y[0]]
xy = odeint(dynamical_system, x0y0, t)
ind = np.argmax(xy[:,0])
gamma = ((xy[ind,0] - fp[0,0])**2 + (xy[ind,1] - fp[0,1])**2)**0.5
x = xy[:,0]
y = xy[:,1]
# x0y0 = [x2[0],y2[0]]
# xy2 = odeint(dynamical_system, x0y0, t)
# x2 = xy2[:,0]
# y2 = xy2[:,1]

print(gamma)

    
fig, ax = plt.subplots(figsize=(7,7))
ax.streamplot(X, Y, dxdt_streamplot, dydt_streamplot, density = 2)
ax.plot(x, y, '-', color='red', linewidth=2)
ax.plot(x2, y2, '-', color='magenta', linewidth=2)
ax.plot(fp[0,0],fp[0,1], '.', color='black', markersize=15, label='Saddle node')
ax.plot(fp[1,0],fp[1,1], '.', color='magenta', markersize=15, label='Unstable spiral')
ax.set_title('$2.4: \mu = 0.07$')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim(xmin,xmax)
ax.set_ylim(xmin,xmax)
ax.set_box_aspect(1) 

plt.legend(loc="upper left")
plt.savefig('24.a_mu0.07.png', bbox_inches='tight')
plt.show()