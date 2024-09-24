import matplotlib.pyplot as plt
import numpy as np

X = np.array([0.7, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
Y = np.array([1.4734, 1.5375, 1.5642, 1.5878,
              1.6084, 1.6270, 1.6436, 1.6586])
Z = np.array([237, 701, 815, 779, 666, 543, 435, 351])

fig, ax = plt.subplots()

x = Y
y = Z

ax.plot(x, y, label='Q($\lambda_{res}$)')
ax.plot(1.5500, 755, **{'color': 'red', 'marker': 'o'}, label='')
ax.set_xlabel('Resonant wavelength, $\lambda_{res}$, um', size=16)
ax.set_ylabel('Quality factor, $\it{Q}$, a.u.', size=16)
ax.set_title('$a_{start}$ = 0.7 um, $a_{end}$ = 0.6 um, $s_{cav}$ = 0.2 um')
ax.grid()
plt.legend()
plt.show()


fig, ax = plt.subplots()

x = X
y = Z

ax.plot(x, y, label='Q(w)')
ax.plot(0.947, 755, **{'color': 'red', 'marker': 'o'}, label='')
ax.set_xlabel('Width of resonator, $\it{w}$, um', size=16)
ax.set_ylabel('Quality factor, $\it{Q}$, a.u.', size=16)
ax.set_title('$a_{start}$ = 0.7 um, $a_{end}$ = 0.6 um, $s_{cav}$ = 0.2 um')
ax.grid()
plt.legend()
plt.show()
