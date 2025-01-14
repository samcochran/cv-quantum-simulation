#code to produce animations

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8-deep')
from matplotlib.animation import FuncAnimation


def animate_evolution(xvec, equilibria, snapshots, K, interval=20, ymax=None):
    """
    Parameters:
        xvec: ndarray containing the grid points for the spatial domain
        equilibria: ndarray containing the equilibrium positions (dimension should be the same as the number of snapshot lists)
        snapshots: list of snapshot lists; each snapshot list contains K ndarrays of wavefunction amplitudes. Each of these arrays should have the same length as xvec.
        K: integer number of time steps, corresponds to the length of each snapshot list
        labels: optional, used for plotting
        interval: optional, corresponds to the delay between animation frames.
    """
    N = len(equilibria)
    fig, ax = plt.subplots()
    lines = []
    for _ in range(N):
        lines.append(ax.plot([], [])[0])
    # if snapshots2 is not None: ln2, = ax.plot([], [], label=label2)
    # plt.legend(loc='upper right')

    if ymax is None: max_val = .05
    else: max_val = ymax
    for snapshot_list in snapshots:
        for snapshot in snapshot_list:
            val = np.max(snapshot)
            if val > max_val:
                max_val = val
    # if snapshots2 is not None:
    #     for snapshot in snapshots2:
    #         val = np.max(snapshot)
    #         if val > max_val:
    #             max_val = val

    def init():
        L = xvec[-1] - max(equilibria)
        ax.set_xlim(-L, L)
        ax.set_ylim(-.05, max_val + .05)
        # if snapshots2 is not None: return ln1, ln2,
        return lines

    def update(frame):
        for line, snapshot_list, i in zip(lines, snapshots, list(range(len(lines)))):
            line.set_data(xvec + equilibria[i], snapshot_list[frame])
        return lines
        # if snapshots2 is not None:
        #     ln2.set_data(xvec - equilibrium_positions[1], snapshots2[frame])
        #     return ln1, ln2,
        # return ln1,

    ani = FuncAnimation(fig, update, frames=list(range(K)),
                        init_func=init, blit=True, interval=interval)
    return ani
