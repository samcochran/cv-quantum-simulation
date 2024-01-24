#code to produce animations

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8-deep')
from matplotlib.animation import FuncAnimation


def animate_evolution(K, xvec, snapshots1, snapshots2=None, label1='First', label2='Second', interval=20):
    fig, ax = plt.subplots()
    lines = []
    ln1, = ax.plot([], [], label=label1)
    if snapshots2 is not None: ln2, = ax.plot([], [], label=label2)
    plt.legend(loc='upper right')

    max_val = .05
    for snapshot in snapshots1:
        val = np.max(snapshot)
        if val > max_val:
            max_val = val
    if snapshots2 is not None:
        for snapshot in snapshots2:
            val = np.max(snapshot)
            if val > max_val:
                max_val = val

    def init():
        L = xvec[-1]
        ax.set_xlim(-L, L)
        ax.set_ylim(-.05, max_val + .05)
        if snapshots2 is not None: return ln1, ln2,
        return ln1,

    def update(frame):
        ln1.set_data(xvec, snapshots1[frame])
        if snapshots2 is not None:
            ln2.set_data(xvec, snapshots2[frame])
            return ln1, ln2,
        return ln1,

    ani = FuncAnimation(fig, update, frames=list(range(K)),
                        init_func=init, blit=True, interval=interval)
    return ani
