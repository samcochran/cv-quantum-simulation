import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp, simpson
from animate import *
from matplotlib import pyplot as plt
import strawberryfields as sf
from strawberryfields.ops import *
from gate_decompositions import *
import warnings
warnings.filterwarnings(action='ignore', module='strawberryfields')


#set global parameters
m = 1
omega = 1/2
potential = lambda x: 1/2*m*omega**2*x**2 #harmonic oscillator
t = np.linspace(0, 1/omega*2*pi, 100) #run for one full period
# t = np.linspace(0, 1/omega*pi, 100) #run for half a period
x_window = [-3, 3]
x_grid = np.linspace(x_window[0], x_window[1], 1000)
window = np.array([x_window, [-.1, np.max(potential(x_grid))]])
np.random.seed(42)
sf.hbar = 1
eng = sf.Engine('gaussian')


#classical version
def func(t, y):
    q_dot = y[1]/m
    # p_dot = -2*(y[0]**2 - const)*2*y[0]
    p_dot = -m*omega**2*y[0]
    return np.array([q_dot, p_dot])
y0_list = np.random.normal([-2, 0], [1/np.sqrt(2)/np.exp(2), 1/np.sqrt(2)/np.exp(2)], (5000, 2))

classical_snapshots = []
for y0 in y0_list:
    sol = solve_ivp(func, [t[0], t[-1]], y0, t_eval=t)
    classical_snapshots.append(sol.y[0, :])
classical_snapshots = np.array(classical_snapshots).T
hist_bins = np.linspace(window[0, 0], window[0, 1], 100)


#quantum version
def prepare_initial_state(n_qumodes, q, x0=None, p0=None, squeeze=None):
    for i in range(n_qumodes - 1):#the -1 is to account for the ancilla
        Vac | q[i] #initial vacuum state
        if squeeze is not None:
            Sgate(squeeze[i]) | q[i] #squeeze (phase angle zero)
        if x0 is not None:
            Xgate(x0[i]) | q[i]
        if p0 is not None:
            Zgate(p0[i]) | q[i]
    Vac | q[n_qumodes - 1]

def step(dt, q):
    # CXgate(-dt*m*omega**2) | (q[0], q[1])
    # CXgate(dt/m) | (q[1], q[0])
    quadratic2(-m*omega**2, dt, 1, 0, q)
    quadratic2(1/m, dt, 0, 1, q)

def run(n_qumodes, n_steps, dt, xvec, pvec, eng, x0=None, p0=None, squeeze=None):
    prog = sf.Program(n_qumodes)
    with prog.context as q:
        prepare_initial_state(n_qumodes, q, x0, p0, squeeze)
        for i in range(n_steps):
            step(dt, q)
    result = eng.run(prog, shots=1)
    return result.state.x_quad_values(0, xvec, pvec)#hardcoded to just be first qumode for now


# setup, run, and create plots
n_qumodes = 3
dt = t[1] - t[0]
L = 6
p_grid = np.copy(x_grid)*10
quantum_snapshots = []
progress = tqdm(total=len(t))
titles = [r'$t = 0$', r'$t = \omega/4$', r'$t = \omega/2$', r'$t = 3\omega/4$', r'$t = \omega$']
label_count = 0
fig, ax = plt.subplots(1, 5, figsize=(20, 4))
for k in range(len(t)):
    results = run(n_qumodes, k, dt, x_grid, p_grid, eng, x0=[-2, 0], p0=[0, 0], squeeze=[2, 2])
    partition = simpson(y=results, x=x_grid)
    quantum_snapshots.append(results)
    progress.update(1)
    if k == 0 or k == len(t)//4 or k == len(t)//2 or k == len(t)*3//4 or t[k] == t[-1]:
        _, _, bar_container = ax[label_count].hist(classical_snapshots[0, :], hist_bins, density=True, color='lightcoral', label='Classical Simulation', zorder=1)
        line = ax[label_count].plot([], [], c='royalblue', label='Quantum Simulation', zorder=2)[0]
        data = classical_snapshots[k, :]
        n, _ = np.histogram(data, hist_bins, density=True)
        for count, rect in zip(n, bar_container.patches):
            rect.set_height(count*window[1, 1]*.2)
        line.set_data(x_grid, quantum_snapshots[k]*window[1, 1]*.2)
        ax[label_count].plot(x_grid, potential(x_grid), c='k', alpha=.5, label='Potential Function', zorder=0)
        if label_count == 0:
            ax[label_count].legend(fontsize=15, loc='upper right')
        ax[label_count].set_xlim(window[0])
        ax[label_count].set_ylim(window[1])
        ax[label_count].set_title(titles[label_count], fontsize=20)
        ax[label_count].set_xticks([])
        ax[label_count].set_yticks([])
        label_count += 1
plt.tight_layout()
plt.savefig(f'./plots/quadratic.png', dpi=400)
plt.close()
progress.close()


# create animations
ani = animate_evolution(window, x_grid, potential, hist_bins, classical_snapshots, quantum_snapshots, title='Quadratic Potential', interval=20)
ani.save('./animations/quadratic_potential.gif', writer='ffmpeg', dpi=400)
ani.save('./animations/quadratic_potential.mp4', writer='ffmpeg', dpi=400)
plt.close()
