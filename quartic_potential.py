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
const = 1
potential = lambda x: (x**2 - const)**2
t = np.linspace(0, 2, 20)
x_window = [-3, 3]
x_grid = np.linspace(x_window[0], x_window[1], 1000)
window = np.array([x_window, [-1, 1 + np.max(potential(x_grid))]])
np.random.seed(42)
sf.hbar = 1
eng = sf.Engine('fock', backend_options={"cutoff_dim": 15})


#classical version
def func(t, y):
    q_dot = y[1]/m
    p_dot = -2*(y[0]**2 - const)*2*y[0]
    return np.array([q_dot, p_dot])
y0_list = np.random.normal([-2, 0], [1/np.sqrt(2)/np.exp(2), 1/np.sqrt(2)/np.exp(2)], (1000, 2))

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
    quadratic2(4, -dt, 1, 0, q)
    quartic([1, 3, 0, 0], 4, dt, 1, 0, None, None, ancilla=2, q=q)
    quadratic2(1/m, -dt, 0, 1, q)

def run(n_qumodes, n_steps, dt, xvec, pvec, eng, x0=None, p0=None, squeeze=None):
    prog = sf.Program(n_qumodes)
    with prog.context as q:
        prepare_initial_state(n_qumodes, q, x0, p0, squeeze)
        for i in range(n_steps):
            step(dt, q)
    result = eng.run(prog, shots=1)
    return result.state.x_quad_values(0, xvec, pvec)#hardcoded to just be first qumode for now

n_qumodes = 3
dt = t[1] - t[0]
L = 6
p_grid = np.copy(x_grid)
quantum_snapshots = []
progress = tqdm(total=len(t))
for k in range(len(t)):
    results = run(n_qumodes, k, dt, x_grid, p_grid, eng, x0=[-2, 0], p0=[0, 0], squeeze=[2, 2])
    partition = simpson(y=results, x=x_grid)
    quantum_snapshots.append(results/partition)#force normalize; this is a bandaid fix, but the fock backend results in amplitude decay
    progress.update(1)
progress.close()


#create animations
ani = animate_evolution(window, x_grid, potential, hist_bins, classical_snapshots, quantum_snapshots, title='Quartic Potential')
ani.save('./animations/quartic_potential.gif', writer='ffmpeg', dpi=400)
ani.save('./animations/quartic_potential.mp4', writer='ffmpeg', dpi=400)
plt.close()
