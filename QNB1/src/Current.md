Theory:

Write down QME explicitly. Lindblad form is not an explicit master equation.
An explicit master equation must look like

\dv{t}p_i = -\sum_j \omega_{ji}p_i + \sum_j \omega_{ij}p_j

where $\omega{ij}$ is the transition rate from j to i.

Simulation:

Plot trace distance between density and gibbs state. Also plot fidelity. Both as
a function of the bath coupling. Plot for several different kappa. Do all of this
for the 2x2 and 4x4 cases. Compare these cases between Redfield and HEOM.
