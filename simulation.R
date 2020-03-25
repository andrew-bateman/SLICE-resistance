# "simulation.R" runs and plots a host/parasite simulation, with default parameter values, according to the model in:
# Bateman, A.W., Peacock, S.J., Krkosek, M., and Lewis, M.A. 2020. Migratory hosts can maintain the high-dose refuge effect in a structured host-parasite system: the case of sea lice and salmon. Evolutionary Applications.
#     Copyright (C) 2020 Andrew Bateman
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


source('setup sim.R')

#run the ODE simulation (from the deSolve package)
out <- ode(State, Time, func = system, Pars,method='lsode',
    events = list(data = Events))

plot.series()

#simulation convert output to a data frame for easy manipulation/viewing
dat = as.data.frame(out)
