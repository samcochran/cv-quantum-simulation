#code to produce animations

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8-deep')
from matplotlib.animation import FuncAnimation


def animate_evolution_particles(x_grid, potential, snapshots, window, interval=60):
    """
    Parameters:
        x_grid: ndarray containing the grid points for the spatial domain
        equilibria: ndarray containing the equilibrium positions (dimension should be the same as the number of snapshot lists)
        snapshots: list of snapshot lists; each snapshot list contains K ndarrays of wavefunction amplitudes. Each of these arrays should have the same length as x_grid.
        window: plotting window, array of the form [[x_min, x_max], [y_min, y_max]]
        interval: optional, corresponds to the delay between animation frames.
    """
    N = snapshots.shape[1]
    K = snapshots.shape[0]
    fig, ax = plt.subplots()
    lines = []
    ax.plot(x_grid, potential(x_grid))
    for _ in range(N):
        lines.append(ax.plot([], [], c='k', marker='.', linestyle=None)[0])

    def init():
        ax.set_xlim(window[0, 0], window[0, 1])
        ax.set_ylim(window[1, 0], window[1, 1])
        return lines

    def update(frame):
        for line, i in zip(lines, list(range(N))):
            line.set_data([snapshots[frame, i]], [potential(snapshots[frame, i])])
        return lines

    ani = FuncAnimation(fig, update, frames=list(range(K)),
                        init_func=init, blit=True, interval=interval)
    return ani
