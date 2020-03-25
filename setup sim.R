# "setup sim.R" sets up a host/parasite simulation, with default parameter values, according to the model in:
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

library(deSolve)

#default parameter values - see main text for details
farm.hosts = 6e6 # approximately six million fish, given 1.5 year grow-out time, 4.9kg per fish - Frazer et al. 2012
beta = 0.22*1.77e-7
m_out.0 = 1/0.25; m_in.0 = 1/1.25; sigma.0 = 4    #salmon migration (juveniles migrate to sea in about 3 months, adults return after 15 months, spawners spend ~3 months in nearshore)
a_0.0 = 5e7; b.0 = 1e5 
mu_j.0 = 19*m_out.0 
alpha.0 = 7.3; delta.0 = 3                #juvenile-host mortality and density-dependent louse mortality (Peacock et al. 2014 - 0.02 hosts/parasite/day)
mu_s.0 = 6.08; mu_r.0 = 6.38              #susceptible and resistant background louse mortality rates (susceptiple 60 d lifespan)
h.0 = 0.67                                #harverst rate (1.5 y schedule)
gamma_T.0 = 25
epsilon_r.0 = 0.05       #treatment level, susceptibility of resistant lice, 
lambda_s.0 = 2.32e3; lambda_r.0 = 2.32e3; c.0 = 73    #larval production rate for lice, mortality rate for larvae (survive about 5 d)
wild.factor = 20 #size of wild (oceanic, unexposed) salmonid population, relative to farm population
wild.transmission.factor = 1/10

print.error = TRUE

#set up the system of differential equations
system <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        #domestic
        dF_D = 0
        dL_sD = lambda_s/c*(beta_DD*L_sD + beta_DE*L_sU*m_in/sigma*(J_E*m_out/m_in)/(F_U + J_E*m_out/m_in))*F_D - (mu_s+h)*L_sD - gamma_T*(L_sD + L_rD)/F_D*L_sD
        dL_rD = lambda_r/c*(beta_DD*L_rD + beta_DE*L_rU*m_in/sigma*(J_E*m_out/m_in)/(F_U + J_E*m_out/m_in))*F_D - (mu_r+h)*L_rD - epsilon_r*gamma_T*(L_sD + L_rD)/F_D*L_rD
        
        #exposed
        dJ_E = a_0*J_E/(b + J_E) - (mu_J + m_out)*J_E - alpha*(L_sE + L_rE)
        dL_sE = lambda_s/c*(beta_ED*L_sD)*J_E - (mu_r + mu_J + m_out)*L_sE - alpha*L_sE*(1 + (L_sE + L_rE)/J_E)
        dL_rE = lambda_r/c*(beta_ED*L_rD)*J_E - (mu_s + mu_J + m_out)*L_rE - alpha*L_rE*(1 + (L_sE + L_rE)/J_E)
      
        #unexposed
        dF_U = 0
        dL_sU = lambda_s/c*beta_UU*L_sU*(F_U + J_E*m_out/m_in) + m_out*L_sE - (mu_s + delta + m_in)*L_sU - delta*L_sU*(L_sU + L_rU)/(F_U + J_E*m_out/m_in)
        dL_rU = lambda_r/c*beta_UU*L_rU*(F_U + J_E*m_out/m_in) + m_out*L_rE - (mu_r + delta + m_in)*L_rU - delta*L_rU*(L_sU + L_rU)/(F_U + J_E*m_out/m_in)
        
        derivs = c(dF_D, dL_sD, dL_rD,
            dJ_E, dL_sE, dL_rE, 
            dF_U, dL_sU,dL_rU)
        
        derivs[!is.finite(derivs)] = 0
        
        return(list(derivs))
    })
}

n = 5000 # number of simulations
Time = seq(0, 5000, length.out = n+1)

#starting population values (to be simulated to equilibrium before inoculation with resistant lice
State = c(F_D = farm.hosts, L_sD = 0, L_rD = 0,
    J_E = farm.hosts*0.5, L_sE = 0, L_rE = 0, 
    F_U = farm.hosts*wild.factor, L_sU = farm.hosts*2, L_rU = 0) #compare F_U=0 to F_U=20

#set parameters to default values
Pars = c(
    a_0 = a_0.0, b = b.0, mu_J = mu_j.0,        
    m_out = m_out.0, m_in = m_in.0, sigma = sigma.0,    
    alpha = alpha.0, delta = delta.0,            
    mu_s = mu_s.0, mu_r = mu_r.0,              
    h = h.0,                            
    gamma_T = gamma_T.0, epsilon_r = epsilon_r.0,      
    lambda_s = lambda_s.0, lambda_r = lambda_r.0, c = c.0,  
    beta_DD = beta, beta_ED = beta/2, beta_UU = beta*wild.transmission.factor, beta_DE = beta/2)

#inoculate system with one resistant louse at the halfway mark
Events = data.frame(var = c("L_rD"),time = c(Time[n/2+1]),value = 1#c(0.01)
                    ,method = c("add"))



#plotting function
plot.series = function(times=range(dat$time),louse.range.scale = 1){
    p0 = par()
    par(mfrow = c(2,2), mar=c(1.5,3,3,0.5))
    dat = as.data.frame(out)
    
    plot(dat$time, dat$F_D, xlim=times, ylim=c(0,1.2*max(dat$F_D)), type="l", col='red', 
        main = paste('wild population:', max(dat$F_U)), axes=FALSE)
    box(); axis(2); title(sub='t ------>', line=0.5)
    lines(dat$time, dat$J_E, col='blue')
    lines(dat$time, dat$F_U)
    text(n,dat$J_E[n]+0.15*max(dat$F_D),'exposed hosts', col='blue', adj=1)
    text(n,dat$F_D[n]+0.15*max(dat$F_D),'farm hosts', col='red', adj=1)
    
    max.y = max(c(dat$L_sD, dat$L_rD),na.rm=TRUE)
    plot(dat$time, dat$L_sD, xlim=times, 
        ylim=louse.range.scale*c(0,1.2*max.y), 
        type="l", col='red', main = 'lice on domestic hosts', axes=FALSE)
    box(); axis(2); title(sub='t ------>', line=0.5)
    lines(dat$time, dat$L_rD, lwd=2, col='red')
    text(n,dat$L_sD[n]+0.15*max.y,'susceptible', col='red', adj=1)
    text(n,dat$L_rD[n]+0.15*max.y,'resistant', col='red', adj=1)
    
    max.y = max(c(dat$L_sE, dat$L_rE),na.rm=TRUE)
    plot(dat$time, dat$L_sE, xlim=times, 
        ylim=louse.range.scale*c(0,1.2*max.y), 
        type="l", col='blue', main = 'lice on exposed juveniles', axes=FALSE)
    box(); axis(2); title(sub='t ------>', line=0.5)
    lines(dat$time, dat$L_rE, lwd=2, col='blue')
    text(n,dat$L_sE[n]+0.15*max.y,'susceptible', col='blue', adj=1)
    text(n,dat$L_rE[n]+0.15*max.y,'resistant', col='blue', adj=1)
    
    max.y = max(c(dat$L_sU, dat$L_rU),na.rm=TRUE)
    plot(dat$time, dat$L_sU, xlim=times, 
        ylim=louse.range.scale*c(0,1.2*max.y), 
        type="l", col='black', main = 'lice in unexposed environment', axes=FALSE)
    box(); axis(2); title(sub='t ------>', line=0.5)
    lines(dat$time, dat$L_rU, lwd=2, col='black')
    text(n,dat$L_sU[n]+0.15*max.y,'susceptible', adj=1)
    text(n,dat$L_rU[n]+0.15*max.y,'resistant', adj=1)
    
    par = p0
}
