#cv quantum simulation of schrodinger evolution of n particles in 1 dimension
#planned extension to 3 dimensions and liouvillian simulation next

import numpy as np
from numpy import pi, sqrt, sin, cos, exp
import strawberryfields as sf
from strawberryfields.ops import *
# import warnings
# warnings.filterwarnings(action='ignore', module='strawberryfields')
from tqdm import tqdm
np.random.seed(42)
sf.hbar = 1


class quadratic_simulator:

    def __init__(self, F, masses=None):
        """
        Parameters:
            F: nxn array of force constants, n is number of particles in the system
        """
        self.F = F
        self.n = F.shape[0]
        if masses is not None:
            self.masses = masses
        else:
            self.masses = np.ones(self.n)

    def prepare_initial_state(self, q, x0=None, p0=None, squeeze=None):
        for i in range(self.n):
            Vac | q[i] #initial vacuum state
            if x0 is not None:
                Xgate(x0[i]) | q[i] #displace q by x0
            if p0 is not None:
                Zgate(p0[i]) | q[i] #displace p by p0
            if squeeze is not None:
                Sgate(squeeze[i]) | q[i] #squeeze (phase angle zero)

    def square_term(self, q, i, dt):
        """
        Quantum circuit for the term e^(-1j t F_ii x_i^2)
        Parameters:
            i: index of the square term being multiplied in the potential, corresponding to the location of the correct force constant along the diagonal of F
        """
        Pgate(-2*dt*self.F[i, i]) | q[i]

    def cross_term(self, q, i, j, dt):
        """
        Quantum circuit for the term e^(-1j t F_ij x_i x_j)
        Parameters:
            i, j: indices of the terms being multiplied in the potential, corresponding to the location of the correct force constant in F
        """
        CZgate(-dt*self.F[i, j]) | (q[i], q[j])

    def kinetic_step(self, q, dt):
        for i in range(self.n):
            Fouriergate() | q[i]
            Pgate(-dt/self.masses[i]) | q[i]
            Fouriergate().H | q[i]

    def potential_step(self, q, dt):
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    self.square_term(q, i, dt)
                else:
                    self.cross_term(q, i, j, dt)

    def run(self, K, dt, xvec, pvec, eng, x0=None, p0=None, squeeze=None):
        prog = sf.Program(self.n)
        with prog.context as q:
            self.prepare_initial_state(q, x0, p0, squeeze)
            for _ in range(K):
                self.kinetic_step(q, dt)
                self.potential_step(q, dt)

        result = eng.run(prog, shots=1)
        return [result.state.x_quad_values(i, xvec, pvec) for i in range(self.n)]
