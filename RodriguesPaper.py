# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt 


# Define your ODE Model 
# y are state variables (dependent on time and parameters)
# t is for time
# k is for parameters (rates, proportions, etc.)
def model(y,t,k):  
    # State Variables
    N = y[0]
    I = y[1]
    Q = y[2]


    # Parameters (changing variables)
    mu = k[0]
    r = k[1]
    K = k[2]
    c1 = k[3]
    a = k[4]
    lamb = k[5]
    s0 = k[6]
    d = k[7]
    rho = k[8]
    gamma = k[9]
    c2 = k[10]
    delta = k[11]
    b = k[12]
    sc = k[13]
    
    mu,r,N,K,c1,I,a,Q,lamb,s0,d,rho,gamma,c2,delta,b

    St = sc
    #St = 0
    qT = 0


    Dosage = (600*1.8)/(1 / 8) 
    tau = 1/8
    T = 21
    
    if (t>=0) and (t<tau): 
            qT = Dosage
    elif (t>=T) and t<(T+tau):  
            qT = Dosage
    
    elif (t>=2 * T) and t<((2*T)+tau):
            qT = Dosage	
    
    elif t>=(3*T) and t<((3*T)+tau): 
            qT = Dosage
    
    elif t>= (4*T) and t<((4*T)+tau): 
            qT = Dosage
    
    elif t>= (5*T) and t<((5*T)+tau): 
        	  qT = Dosage
    
    else:
            qT=0
    
    print([t, qT])


    # Differential Equations
    dNdt = (r*N)*(1-(N/K))-(c1*N*I)-((mu*N*Q)/(a+Q))
    dIdt = St+s0-(d*I)+((rho*N*I)/(gamma+N))-(c2*N*I)-((delta*I*Q)/(b+Q))
    dQdt = qT-(lamb*Q)
    
    
    # Compile and return    
    dydt = [dNdt, dIdt, dQdt]
    return dydt

# Initial Conditions
N0 = 2e10
I0 = 5e7
Q0 = 0



y0 = [N0, I0, Q0]



# Parameter Values
mu = 8
r = 0.01
K = 1e12
c1 = 5e-11
I = 0.01
a = 2e3
Q = 5 
lamb = 4.16
s0 = 3e5
d = 0.001
rho = 1e-12
gamma = 100
c2 = 1e-13
delta = 10000
b = 5e6
sc = 0

p = [mu,r,K,c1,a,lamb,s0,d,rho,gamma,c2,delta,b,sc]

# Setup Time
tB = 0  #Beginning Time
tF = 2000 # Final Time
nTime = 2000 # time points

# Create a vector of time points from beginning time to end time that includes nTime number of points
tspan = np.linspace(tB,tF,nTime)

ode_out = odeint(model, y0,tspan, args=(p,), rtol=1e-9, atol=1e-9, hmax=1/10)

N = ode_out[:,0]
I = ode_out[:,1]
Q = ode_out[:,2]



# Setup Time
tB = 0  #Beginning Time
tF = 100 # Final Time
nTime = 100 # time points



# Plot Results
plt.semilogy(tspan, I, 'r:', label="Infectious")  # Note you need the labels for the plot legend
plt.ylabel('Infected Cells (i)') # Vertical axis label
plt.xlabel('Time (days)') # Horizontal Axis label
plt.legend(loc='best') # Show legend and put it in the "best" location
plt.ylim( (10**0,10**16) )
plt.show() # Show PLot

