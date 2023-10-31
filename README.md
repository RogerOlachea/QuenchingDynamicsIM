# QuenchingDynamicsIM
Zero temperature quenching dynamics of Ising Model under external magnetic field.

In this repository are the fortran codes for the quenching dynamics of the ising model under external magnetic field, the code uses dirty random number generator from numerical recipes and glauber dynamics algorithm for the flip acceptance ratio, some measures are made like magnetization, persistence, autocorrelation and domain walls number.
In the "Data" folder are the measurement data for lattice sizes 100, 300, 500 and 800.

IM stands for IM interaction without external field (h = 0)
IMUF for IM interaction with uniform field (h = +1)
IMRUF for IM interaction with random uniformly distributed between -1 and 1 field (h = [-1,1])
IMRBF for IM interaction with random binary field (h = {-1,1})
IMRNF for IM interaction with random normaly distributed betweeen -3.7 and 3.7 field (h = [-3.7,3.7])

