# CV-Quantum Simulation of Classical Dynamics

The code in this repository can be used to simulate classical Hamiltonian systems using the Koopman von Neumann formalism, as described in our paper [An application of continuous-variable gate synthesis to quantum simulation of classical dynamics](https://arxiv.org/abs/2407.08006). Gate decompositions are implemented for operations up to fourth order in the quadrature operators, though technological limitations lead to challenges in implementing the higher order terms. 

Below are animations produced using our code for two different potential functions. For the quadratic potential, which corresponds to a harmonic oscillator, the quantum simulation matches well with the classical molecular dynamics simulation.

<img src="/animations/quadratic_potential.gif " width="600"/>

The inaccuracies in the simulation of the quartic potential arise from the truncation of the Fock basis in the simulation backend. Decompositions for cubic and quartic terms rely heavily on the non-Gaussian cubic phase gate, necessitating the use of a Fock backend which truncates the Hilbert space. We have found that the cubic phase gate is particularly prone to these cutoff errors. While the simulation of the quartic potential is at least directionally correct, it is clearly unable to accurately describe the dynamics. Future work will involve exploring alternative simulation platforms to try to address this issue.

<img src="/animations/quartic_potential.gif " width="600"/>

