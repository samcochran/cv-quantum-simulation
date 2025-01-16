import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint, quad
from animate_distribution import animate_evolution
from animate_particles import animate_evolution_particles

#potential function is (x_1**2 + x_2**2 - 1)**2
# H = 1/(2*m)*p**2 + (x**2 - 1)**2

m = 1
potential = lambda x: (x**2 - 10)**2
x_window = [-5, 5]
x_grid = np.linspace(x_window[0], x_window[1], 1000)
window = np.array([x_window, [-1, 1 + np.max(potential(x_grid))]])

def func(y, t):
    q_dot = y[1]/m
    p_dot = -2*(y[0]**2 - 1)*2*y[0]
    return np.array([q_dot, p_dot])
t = np.linspace(0, 5, 250)
squeeze_factor = 10
y0_list = np.random.normal([-4, -5], [1/np.sqrt(2)/squeeze_factor, 1/np.sqrt(2)*squeeze_factor], (1000, 2))
# plt.plot(x_grid, potential(x_grid))
# plt.scatter(y0_list[:, 0], potential(y0_list[:, 0]))
# plt.hist(y0_list[:, 0], range=window[0], bins=100, density=True)
# plt.show()

snapshots = []
for y0 in y0_list:
    sol = odeint(func, y0, t)
    snapshots.append(sol[:, 0])
snapshots = np.array(snapshots).T

hist_bins = np.linspace(window[0, 0], window[0, 1], 100)
ani = animate_evolution(x_grid, potential, hist_bins, snapshots, window)
ani.save('./animations/quartic_potential_hist.mp4', writer='ffmpeg')
plt.close()

ani = animate_evolution_particles(x_grid, potential, snapshots, window)
ani.save('./animations/quartic_potential_particles.mp4', writer='ffmpeg')
plt.close()
