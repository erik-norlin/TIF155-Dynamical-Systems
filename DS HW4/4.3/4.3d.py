# Exercise 4.3d

import numpy as np
import matplotlib.pyplot as plt
import sys

title = '4.3d'

a = 1.4
b = 0.3

cut_tail = 100

T = 10**4
dt = 5*10**-3
t = np.arange(0,T,dt)

x = np.zeros_like(t,dtype=float)
y = x.copy()

x[0] = (np.random.uniform()-0.5)
y[0] = (np.random.uniform()-0.5)

for t in range(len(x)-1):
    x[t+1] = y[t]+1-a*x[t]**2
    y[t+1] = b*x[t]

x = x[cut_tail:]
y = y[cut_tail:]

q_vals = np.linspace(0,4,10)
Dq_list = []
epsilon_range = np.linspace(10**-3, 2*10**-2, 10)
fig, ax = plt.subplots(figsize=(7,7))
color = ['tab:blue', 'tab:green', 'tab:red']

for q_i in range(len(q_vals)):
    
    q = q_vals[q_i]
    Iq = np.zeros_like(epsilon_range, dtype=float)

    for epsilon_i in range(len(epsilon_range)):
        
        epsilon = epsilon_range[epsilon_i]
        N_points = len(x)

        xmax = 1.3
        xmin = -xmax
        ymax = 0.4
        ymin = -ymax

        x_bins = np.linspace(xmin, xmax, int((xmax-xmin)/epsilon))
        y_bins = np.linspace(ymin, ymax, int((ymax-ymin)/epsilon))

        plt.figure()
        histogram = plt.hist2d(x, y, bins=[x_bins, y_bins])
        boxes = histogram[0].copy()
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
        
        for i in range(len(boxes[:,0])):
            for j in range(len(boxes[0,:])):
                if boxes[i,j] != 0:
                    if q == 1:
                        N_k = boxes[i,j]
                        Iq[epsilon_i] += ((N_k/N_points)*np.log(1/(N_k/N_points)))
                    else:
                        N_k = boxes[i,j]
                        Iq[epsilon_i] += (N_k/N_points)**q

    x_axis = np.log(1/epsilon_range)

    if q == 1:
        y_axis = Iq
    else:
        y_axis = np.log(Iq)/(1-q)

    coef = np.polyfit(x_axis, y_axis, 1)
    Dq = coef[0]
    Dq_list.append(Dq)
    print('Dq=',Dq,'    q=',q)

ax.plot(q_vals, Dq_list, 'o--', color='black')
ax.set_title(title)
ax.set_box_aspect(1)
ax.set_xlabel(r'$q$')
ax.set_ylabel(r'$D_{q}$')

fig.savefig('Dynamical Systems/DS HW4/4.3/'+title+'.png')
plt.show()