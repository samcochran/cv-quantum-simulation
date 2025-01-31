#code to produce animations

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8-deep')
from matplotlib.animation import FuncAnimation

def animate_evolution(window, x_grid, potential, hist_bins, classical_snapshots, quantum_snapshots=None, interval=60, title=None):
    """
    Parameters:
        x_grid: ndarray containing the grid points for the spatial domain
        equilibria: ndarray containing the equilibrium positions (dimension should be the same as the number of snapshot lists)
        snapshots: list of snapshot lists; each snapshot list contains K ndarrays of wavefunction amplitudes. Each of these arrays should have the same length as x_grid.
        window: plotting window, array of the form [[x_min, x_max], [y_min, y_max]]
        interval: optional, corresponds to the delay between animation frames.
    """
    N = len(classical_snapshots)
    K = classical_snapshots.shape[0]
    n, _ = np.histogram([], hist_bins)
    fig, ax = plt.subplots()
    if title is not None: fig.suptitle(title, fontsize=14)
    ax.plot(x_grid, potential(x_grid), c='k', alpha=.5, label='Potential Function', zorder=0)
    _, _, bar_container = ax.hist(classical_snapshots[0, :], hist_bins, density=True, color='lightcoral', label='Classical Simulation', zorder=1)
    line = ax.plot([], [], c='royalblue', label='Quantum Simulation', zorder=2)[0]
    ax.legend()

    def init():
        ax.set_xlim(window[0, 0], window[0, 1])
        ax.set_ylim(window[1, 0], window[1, 1])
        return bar_container

    def update(frame):
        data = classical_snapshots[frame, :]
        n, _ = np.histogram(data, hist_bins, density=True)
        for count, rect in zip(n, bar_container.patches):
            rect.set_height(count*window[1, 1]*.1)
            if quantum_snapshots is not None:
                line.set_data(x_grid, quantum_snapshots[frame]*window[1, 1]*.1)
        return bar_container.patches

    ani = FuncAnimation(fig, update, frames=list(range(K)),
                        init_func=init, blit=True, interval=interval)
    return ani
