# CV-Quantum Simulation of Classical Dynamics

The code in this repository implements a continuous-variable quantum algorithm to simulate classical Hamiltonian systems using the Koopman von Neumann formalism, as described in our paper [An application of continuous-variable gate synthesis to quantum simulation of classical dynamics](https://doi.org/10.1116/5.0234007). The algorithm is simulated using Strawberry Fields. Gate decompositions are implemented for operations up to fourth order in the quadrature operators, though technological limitations lead to challenges in implementing the higher order terms. The animations below were produced by applying our code to two different potential functions. 

For a quadratic potential corresponding to a simple harmonic oscillator, the quantum simulation matches well with the classical molecular dynamics simulation.

<img src="/animations/quadratic_potential.gif " width="750"/>

Simulation of higher-degree potentials leads to large numerical inaccuracies, which arise from the truncation of the Fock basis in the simulation backend. In particular, decompositions for cubic and quartic terms rely heavily on the non-Gaussian cubic phase gate; as noted in the [Strawberry Fields documentation](https://strawberryfields.ai/photonics/conventions/gates.html), the existing implementation of the cubic phase gate is known to suffer heavily from these inaccuracies. While we observe convergence as the number of Fock basis elements increases, the convergence is slow, especially when other gates are involved alongside the cubic phase, making simulations of higher-degree systems computationally impractical with currently available resources. The inaccuracies become more extreme whenever two-qumode gates are needed on modes where cubic phases are also applied.

We analyzed the cubic phase error by applying a single cubic phase gate to the vacuum state, followed by a two-qumode identity operation on a second mode. The results are summarized in the following observations:
1. Mean Squared Error (MSE): The MSE of the Born rule probability distribution decreases as the cutoff dimension increases, but it remains larger than we need for all tested cutoff dimensions.
2. Density Matrix Trace: The trace of the density matrix is less than one for all tested cutoff dimensions, indicating a nonphysical state.
3. Largest Eigenvalue: The largest eigenvalue of the density matrix is also less than one, suggesting spurious correlations with the second qumode.

![Image](https://github.com/user-attachments/assets/8a0f7019-3580-4830-b058-55134623b3fc)

When multiple cubic phase gates must be applied in sequence, these errors quickly accumulate. For example, when trying to implement a decomposition for the quartic phase gate, which involves multiple cubic phase gates along with two-mode operations with an ancilla qumode, the errors are severe enough to prevent accurate simulation of even a single quartic phase gate.

![Image](https://github.com/user-attachments/assets/4f9cb68d-1e83-49cb-bf1d-4cf80026c679)

When increasing the cutoff dimension past 70 in an effort to increase accuracy, overflow errors prevent the simulation from even running.
