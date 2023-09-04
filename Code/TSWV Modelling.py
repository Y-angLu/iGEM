#Waterloo iGEM 2023

#Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

#Model without time delay or age structures 

#Initial conditions

vh=5 #Birth rate of host
Bv=0.6 #Bite rate of thrip
Thv=0.4 #Infection rate from vector to host
yh=40 #Age at harvest of host
uh=0.015 #Death rate of host 
th=4 #Incubation period for exposed host 
vv=5 #Birth rate of host 
Kv=10000 #Capacity of thrip population
Tvh=0.2 #Infection transmission rate from host to vector 
uv=0.03 #Death rate of thrip adults 
ul=0.4 #Death rate of thrip larvae
yv=10 #Age at maturity of WFT
tv=4 #Incubation period for vector 

def hostpopulation():
    Sh=80 #number of susceptible hosts
    Eh=0 #number of exposed hosts 
    Ih=1 #number of infected hosts
    Iv=0 #number of infective hosts 
    Hh=0 #number of harvest hosts
    Nh=90 #total number of hosts 
    dayplot=[]
    Shplot=[]
    Ehplot=[]
    Ihplot=[]
    Hhplot=[]
    Nhplot=[]
    Ivplot=[]
    for i in range(1,181):
            dayplot.append(i)
            Shplot.append(Sh)
            Ehplot.append(Eh)
            Ihplot.append(Ih)
            Hhplot.append(Hh)
            Nhplot.append(Nh)
            Ivplot.append(Iv)
            dSh = vh - (Bv*Thv*Iv*Sh)/Nh - (1/yh)*Sh - uh*Sh
            dEh = (Bv*Thv*Iv*Sh)/Nh - (1/th)*Eh - (1/yh)*Eh - uh*Eh
            dIh = (1/th)*Eh - (uh*Ih)
            dHh = ((1/yh)*Sh) + ((1/yh)*Eh)
            dNh = (vh) - (uh*Sh) - (uh*Eh) - (uh*Ih)
            dIv = ((1/yv)*Ih) - (uv*Ih)
            Sh += dSh
            Eh += dEh
            Ih += dIh
            Hh += dHh
            Nh += dNh
            Iv += dIv
    plt.plot(dayplot,Shplot)
    plt.plot(dayplot,Ehplot)
    plt.plot(dayplot,Ihplot)
    plt.plot(dayplot,Hhplot)
    #plt.plot(dayplot,Nhplot)
    #plt.plot(dayplot,Ivplot)
    plt.show()

Sl=0 #number of susceptible larva 
El=0 #number of exposed larva
Il=0 #number of infected larva
Sv=0 #number of non-infective adult
Iv=0 #number of infective adult
Nv=0 #total number of WFT

#def WFTpopulation():
#    dSl = vv*(Sv+Iv)*(1-(Nv/Kv)) - (Bv*Tvh*Ih*Sl)/Nh - (1/yv)*Sl - ul*Sl
#    dEl = (Bv*Tvh*Ih*Sl)/Nh - (1/tv)*El - (1/yv)*El - ul*El
#    dIl = (1/tv)*El - (1/yv)*Il - ul*Il
#    dSv = (1/yv)*Sl + (1/yv)*El - uv*Sv
#    dIv = (1/yv)*Il - uv*Iv
#    dNv = vv*(Sv+Iv)*(1-(Nv/Kv)) - uv*Nv

hostpopulation()
