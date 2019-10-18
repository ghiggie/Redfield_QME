In this folder, I will experiment with my code in order to ensure that it returns
consistent results. I will study first a single qubit with Hamiltonian
$H_S = \frac{\omega_0}{2}\sigma_z$, coupled to the bath through the operator
$V_I = \sigma_x$. I will run the code for $100$ time steps, at a bath coupling
strength of 0.1. This isn't a small enough coupling strength for the Redfield
equation to be valid, but I just want to check the results of the program.

I will then study a two qubit system such that the two qubits are not
interacting. As such, we can write the Hamiltonian of the qubits as
$H_S = \mathscr{I}\bigotimes\frac{\omega_0}{2}\sigma_z$. I do not want the left
qubit to interact with the bath, so I will write the coupling operator as
$V_I = \mathscr{I}\bigotimes\sigma_x$. As before, I will run the code for $100$
time steps at a bath coupling strength of 0.1. However, in this case, I will
then use partial trace to extract the dynamics of the second qubit. Because the
first qubit doesn't actually interact with the other objects, this dynamics
should be identical to the dynamics obtained in the first case.

Note: it will be interesting to see what the entropy of the total and partial
systems are in the second case. Theoretically, I expect that the entropy of the
first qubit will be constant and $0$. This will be a good testing ground for my
program.
