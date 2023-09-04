#Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Bacteria: Pseudomonas Aeruginosa -> Logistic Growth Assumed
#Article: https://www.frontiersin.org/articles/10.3389/fsysb.2022.899990/full

#Variables 
    #'H' -> Rod-shaped bacteria with no internalised AB:
    #'I' -> Rod-shaped bacteria with internalised AB:
    #'L' -> Spherical (L-form) cells
    #'P' -> Phagocytes

#Subscripts:
    #'S' -> Antibiotic-susceptible
    #'R' -> Antibiotic-resistant

#Initial Conditions/Values:

#[Hs, Hr, Is, Ir, Ls, Lr, P]
init_values = [10**4, 0, 0, 0, 0, 0, 0]

r = 3 # -> Specific growth rate
c = 0.15 # -> Fitness cost parameter
k = 10**6 # -> Carrying capacity

v = 1.57 # -> Maximum AB internalisation rate 
w = 0.5 # -> Efflux pump upregulator parameter
A50 = 0.47 # -> AB concentration required for half maximal internalisation

y = 1.57 # -> Maximum transition rate
T50 = 0.47 # -> Concentration needed for half-maximal transition
ki = 0.5 # -> Maximum recruitment rate of phagocytes by bacteria

PMAX = 1.8*(10**5) # -> Carrying capacity of phagocytes 
m = 0.017 # -> Death rate of spherical cells 
p = 100 # -> Death rate of rod-shaped bacteria due to AB

oi = 2.4**(10**(-4)) # -> Death rate of bacteria by phagocytes 
u = 6**(10**(-6)) # -> Death rate of phagocytes by bacteria-induced apoptosis 
E = 1.512 # -> Natural clearance rate of phagocytes

a = 15 # -> Extracellular AB decay rate
vl = 1**(10**(-3)) # -> #xtracellular AB decay due to internalisation in rod-shaped cells 
A = 10 # -> Extracellular AB concentration 

#Differential Equations:
def deq(init_values, t):
    Hs, Hr, Is, Ir, Ls, Lr, P = init_values

    dHs_dt = r*Hs*(1-((Hs+Hr)/k)) - ((v*A)/(A+A50))*Hs - ((y*A)/(A+T50))*Hs - oi*Hs*P
    dHr_dt = (1-c)*r*Hr*(1-((Hs+Hr)/k)) - (1-w)*((v*A)/(A+A50))*Hr - ((y*A)/(A+T50))*Hr - oi*Hr*P
    dIs_dt = ((v*A)/(A+A50))*Hs - p*Is - oi*Is*P
    dIr_dt = (1-w)*((v*A)/(A+A50))*Hr - p*Ir - oi*Ir*P
    dLs_dt = ((y*A)/(A+T50))*Hs - m*Ls - oi*Ls*P
    dLr_dt = ((y*A)/(A+T50))*Hr - m*Lr - oi*Lr*P
    dP_dt = ki*(Hs+Is+Hr+Ir)*(1-(P/PMAX)) + ki*(Ls+Lr)*(1-(P/PMAX)) - u*(Hs+Is+Hr+Ir+Ls+Lr)*P - E*P

    return[dHs_dt, dHr_dt, dIs_dt, dIr_dt, dLs_dt, dLr_dt, dP_dt]

#Time
t = np.linspace(0,7.5,1000)

#Solving differential equations
sE = odeint(deq, init_values, t)

# Solutions for plotting
Hs, Hr, Is, Ir, Ls, Lr, P = sE[:, 0], sE[:, 1], sE[:, 2], sE[:, 3], sE[:, 4], sE[:, 5], sE[:,6]

# Plotting
plt.figure(figsize=(8, 6))

plt.plot(t, Hs, label='Hs(t)')
plt.plot(t, Hr, label='Hr(t)')
plt.plot(t, Is, label='Is(t)')
plt.plot(t, Ir, label='Ir(t)')
plt.plot(t, Ls, label='Ls(t)')
plt.plot(t, Lr, label='Lr(t)')

plt.xlabel('Time (days)')
plt.ylabel('Cell Growth (cm-3)')

plt.legend()
plt.grid(True)
plt.show()
