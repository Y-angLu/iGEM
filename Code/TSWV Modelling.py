#Libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint

#Modelling tomato spotted wilt virus -> no time delay or age structures 
#Article: https://www.researchgate.net/publication/271702772_New_mathematical_models_for_vector-borne_disease_transmission_of_tomato_spotted_wilt_virus

#Initial conditions/values

#Host population:
#Sh -> number of susceptible hosts
#Eh -> number of exposed hosts 
#Ih -> number of infected hosts
#Hh -> number of harvest hosts
#Nh -> total number of hosts 

#WFT population:
# Sl -> number of susceptible larva
# El -> number of exposed larva
# Il -> number of infected larva
# Sv -> number of non-infective larva
# Iv -> number of infective adult
# Nv -> total number of WFT

#[Sh, Eh, Ih, Hh, Nh, Sl, El, Il, Sv, Iv, Nv]
init_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

vh = 5 #Birth rate of host
Bv = 0.6 #Bite rate of thrip
Thv = 0.4 #Infection rate from vector to host
yh = 40 #Age at harvest of host
uh = 0.015 #Death rate of host 
th = 4 #Incubation period for exposed host 
vv = 5 #Birth rate of host 
Kv = 10000 #Capacity of thrip population
Tvh = 0.2 #Infection transmission rate from host to vector 
uv = 0.03 #Death rate of thrip adults 
ul = 0.4 #Death rate of thrip larvae
yv = 10 #Age at maturity of WFT
tv = 4 #Incubation period for vector 

def hostpopulation(init_values, t):
    Sh, Eh, Ih, Hh, Nh, Sl, El, Il, Sv, Iv, Nv = init_values

    dSh_dt = vh - (Bv*Thv*Iv*Sh)/Nh - (1/yh)*Sh - uh*Sh
    dEh_dt = (Bv*Thv*Iv*Sh)/Nh - (1/th)*Eh - (1/yh)*Eh - uh*Eh
    dIh_dt = (1/th)*Eh - (uh*Ih)
    dHh_dt = (1/yh)*Sh + ((1/yh)*Eh)
    dNh_dt = (vh) - (uh*Sh) - (uh*Eh) - (uh*Ih)
        
    dSl_dt = vv*(Sv+Iv)*(1-(Nv/Kv)) - (Bv*Tvh*Ih*Sl)/Nh - (1/yv)*Sl - ul*Sl
    dEl_dt = (Bv*Tvh*Ih*Sl)/Nh - (1/tv)*El - (1/yv)*El - ul*El
    dIl_dt = (1/tv)*El - (1/yv)*Il - ul*Il
    dSv_dt = (1/yv)*Sl + (1/yv)*El - uv*Sv
    dIv_dt = (1/yv)*Il - uv*Iv
    dNv_dt = vv*(Sv+Iv)*(1-(Nv/Kv)) - uv*Nv

    return[dSh_dt, dEh_dt, dIh_dt, dHh_dt, dNh_dt, dSl_dt, dEl_dt, dIl_dt, dSv_dt, dIv_dt, dNv_dt]

t = np.linspace(0, 200, 1000)  
 
sE = odeint(hostpopulation, init_values, t)

Sh, Eh, Ih, Hh, Nh, Sl, El, Il, Sv, Iv, Nv = sE[:, 0], sE[:, 1], sE[:, 2], sE[:, 3], sE[:, 4], sE[:, 5], sE[:, 6], sE[:, 7], sE[:, 8], sE[:, 9], sE[:, 10]

# Plotting
plt.figure(figsize=(8, 6))

#Population dynamics of hosts

# plt.plot(t, Sh, label='Sh(t)')
# plt.plot(t, Eh, label='Eh(t)')
# plt.plot(t, Ih, label='Ih(t)')
# plt.plot(t, Hh, label='Hh(t)')

#Population dynamics of WFT

plt.plot(t, Nh, label='Sl(t)')
plt.plot(t, Sh, label='El(t)')
plt.plot(t, Eh, label='Il(t)')
plt.plot(t, Ih, label='Sv(t)')
plt.plot(t, Hh, label='Iv(t)')

plt.xlabel('Time (days)')
plt.ylabel('Growth')

plt.legend()
plt.grid(True)
plt.show()

