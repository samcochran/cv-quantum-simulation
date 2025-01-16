#code to produce animations

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8-deep')
from matplotlib.animation import FuncAnimation

def animate_evolution(x_grid, potential, hist_bins, snapshots, window, interval=60):
    """
    Parameters:
        x_grid: ndarray containing the grid points for the spatial domain
        equilibria: ndarray containing the equilibrium positions (dimension should be the same as the number of snapshot lists)
        snapshots: list of snapshot lists; each snapshot list contains K ndarrays of wavefunction amplitudes. Each of these arrays should have the same length as x_grid.
        window: plotting window, array of the form [[x_min, x_max], [y_min, y_max]]
        interval: optional, corresponds to the delay between animation frames.
    """
    N = len(snapshots)
    K = snapshots.shape[0]
    n, _ = np.histogram([], hist_bins)
    fig, ax = plt.subplots()
    _, _, bar_container = ax.hist(snapshots[0, :], hist_bins, density=True, color='k', alpha=0.75)
    ax.plot(x_grid, potential(x_grid))

    def init():
        ax.set_xlim(window[0, 0], window[0, 1])
        ax.set_ylim(window[1, 0], window[1, 1])
        return bar_container

    def update(frame):
        data = snapshots[frame, :]
        n, _ = np.histogram(data, hist_bins)
        for count, rect in zip(n, bar_container.patches):
            rect.set_height(count)
        return bar_container.patches

    ani = FuncAnimation(fig, update, frames=list(range(K)),
                        init_func=init, blit=True, interval=interval)
    return ani
