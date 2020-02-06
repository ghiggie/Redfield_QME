Assigned Work
-------------

Theory:

Write down QME explicitly. Lindblad form is not an explicit master equation.
An explicit master equation must look like

\dv{t}p_i = -\sum_j \omega_{ji}p_i + \sum_j \omega_{ij}p_j

where $\omega{ij}$ is the transition rate from j to i.

Simulation:

Plot trace distance between density and gibbs state. Also plot fidelity. Both as
a function of the bath coupling. Plot for several different kappa. Do all of this
for the 2x2 and 4x4 cases. Compare these cases between Redfield and HEOM.

My Work
-------

I will consider first a single qubit coupled to the bath. The Hamiltonian of the
system is

H = \frac{\omega_0}{2}\sigma_z

and the coupling operator to the bath is

V = \lambda \sigma_y

Cases to consider:
1. First, consider the case in which T > \omega_0. Then plot over \lambda = 0.01
   to 0.10. Use various initial states.
2. Consider the case in which T < \omega_0. Do the same as above.
