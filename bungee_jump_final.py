from pylab import *
import numpy as np

dt = 0.01
simulationTime = 60.0
n = int(round(simulationTime / dt))
m = 70.0  # kg
k = 150.0  # N/m
d = 20.0  # m
g = 9.81  # m/s^2
h = 50.0  # m height of bridge
# we find Cv experimentally
Cv = 12.72  # kg/s

t = zeros(n, float)
v = zeros(n, float)
x = zeros(n, float)
Etotal = m * g * h
Ug = zeros(n, float)
Us = zeros(n, float)
Ek = zeros(n, float)

x[0] = 0
v[0] = 0
Ug[0] = Etotal
Us[0] = 0.0
Ek[0] = 0.0

for i in range(0, n - 1):
    if x[i] > d:
        F = -k * (x[i] - d) - Cv * v[i] + m * g
    else:
        F = m * g

    a = F / m

    v[i + 1] = v[i] + a * dt
    x[i + 1] = x[i] + v[i + 1] * dt
    Ug[i + 1] = m * g * (h - x[i + 1])
    Ek[i + 1] = 0.5 * m * v[i + 1] * v[i + 1]

    if (x[i + 1] > d):
        Us[i + 1] = 0.5 * k * (x[i + 1] - d) * (x[i + 1] - d)
    else:
        Us[i + 1] = 0.0

    t[i + 1] = t[i] + dt

xmax = np.max(x)
print("Max distance from the origin: ", xmax)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
fig.canvas.set_window_title('Bungee Jump')

ax1.plot(t, v, 'r')
ax1.set(ylabel='velocity [m/s]')

ldistanceFromOrigin, = ax2.plot(t, x, 'c')
ax2.set(ylabel='distance [m]')


l1, = ax3.plot(t, Ug, 'y')
l2, = ax3.plot(t, Ek, 'g')
l3, = ax3.plot(t, Us, 'b')
l4, = ax4.plot(t, Ug+Ek+Us,'r')
l5, = ax4.plot(t, Ug+Ek, 'b')

ax3.set(ylabel='energy [J]')
ax3.legend((l1, l2, l3), ('gravitational potential energy', 'kinetic energy', 'elastic potential energy'), loc='upper right', prop={'size': 7})

ax4.set(xlabel='time [s]', ylabel='energy [J]')
ax4.legend((l4, l5), ('GPE + KE + EPE', 'GPE + KE'), loc='upper right', prop={'size': 8})

plt.show()
