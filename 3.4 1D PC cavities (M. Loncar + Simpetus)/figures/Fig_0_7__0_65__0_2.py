import matplotlib.pyplot as plt
import numpy as np

X = np.arange(0.7, 1.2, 0.1)
Y = np.array([1.5105, 1.5460, 1.5777, 1.6051, 1.6291])
Z = np.array([231, 366, 458, 492, 479])

# Plot the surface
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(X, Y, zs=Z, zdir='Z', label='Curve')

# ax.legend()
# ax.set_xlabel('Width of resonator, w, um')
# ax.set_ylabel('Resonant wave length, $\lambda_{res}$, um')
# ax.set_zlabel('Quality factor, Q, a.u.')

# plt.show()


fig, ax = plt.subplots()

x = Y
y = Z

ax.plot(x, y, label='Q($\lambda_{res}$)')
ax.plot(1.55, 378, **{'color': 'red', 'marker': 'o'}, label='')
ax.set_xlabel('Resonant wavelength, $\lambda_{res}$, um', size=16)
ax.set_ylabel('Quality factor, $\it{Q}$, a.u.', size=16)
ax.set_title('$a_{start}$ = 0.7 um, $a_{end}$ = 0.65 um, $s_{cav}$ = 0.2 um')
ax.grid()
plt.legend()
plt.show()


fig, ax = plt.subplots()

x = X
y = Z

ax.plot(x, y, label='Q(w)')
ax.plot(0.813, 378, **{'color': 'red', 'marker': 'o'}, label='')
ax.set_xlabel('Width of resonator, $\it{w}$, um', size=16)
ax.set_ylabel('Quality factor, $\it{Q}$, a.u.', size=16)
ax.set_title('$a_{start}$ = 0.7 um, $a_{end}$ = 0.65 um, $s_{cav}$ = 0.2 um')
ax.grid()
plt.legend()
plt.show()
