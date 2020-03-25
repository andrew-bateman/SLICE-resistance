To simulate the host/parasite system from:
    Bateman, A.W., Peacock, S.J., Krkosek, M., and Lewis, M.A. 2020. Migratory hosts 
    can maintain the high-dose refuge effect in a structured host-parasite system: 
    the case of sea lice and salmon. Evolutionary Applications,
at default parameter values, run 'simulation.R' in R.
(Package deSolve needs to be installed first.)

To explore effect of parameters and starting conditions on the system, 
adjust the values of the named vectors "Pars" and "State," respectively.

For extremely slow rates of resistance spread, increase 'n' (maximum time steps) to ensure final results approximate equilibrium.