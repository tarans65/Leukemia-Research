# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 11:39:19 2023

@author: ramka
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
    N1 = y[0]
    N2 = y[1]
    I = y[2]
    Q = y [3]


    # Parameters (changing variables)
    r1 = k[0]
    k1 = k[1]
    alpha12 = k[2]
    c1 = k[3]
    mu = k[4]
    a = k[5]
    r2 = k[6]
    alpha21 = k[7]
    v = k[8]
    b = k[9]
    s = k[10]
    m = k[11]
    rho = k[12]
    gamma = k[13]
    delta = k[14]
    c2 = k[15]
    lamb = k[16]
    c = k[17]
    

# Define q(t)

    if (t>=0) and (t<1/8): # First infusion 
        qT = 8400
        #print(t)
    elif (t>=21) and (t<(21+(1/8))):  # second infusion
        qT = 8400
        #print(t)
    elif t>=42 and t<(42+1/8): # third infusion
        qT = 8400	

    elif t>= 63 and t<(63+1/8): #fourth infusion
        qT = 8400

    elif t>= 80 and t<(80+1/8): #fifth infusion
        qT = 8400
    elif t>= 105 and t<(105+1/8): #sixth infusion
    	qT = 8400
    else:
        qT=0

    #print([t, qT])


    # Differential Equations
    dN1dt = 0
    dN2dt = r2-((v*N2*Q)/(b+Q))
    dIdt = s-(m*I)-((delta*I*Q)/(c+Q))
    dQdt = Q-(lamb*Q)
        
    # Compile and return    
    dydt = [dN1dt, dN2dt, dIdt, dQdt]
    return dydt

# Initial Conditions
N10 = 0
N20 = 1e12
I0 = 1e7
Q0 = 0


y0 = [N10, N20, I0, Q0]


# Parameter Values
r1 = 1e-2
k1 = 1e12
alpha12 = 9e-5
c1 = 5e-11
mu = 8
a = 2e3
r2 = 1e7
alpha21 = 9e-16
v = 8
b = 5e6
s = 7e5
m = 1e-3
rho = 1e-12
gamma = 1e2
delta = 1e4
c2 = 1e-13
lamb = 4.16
c = 5e6


p = [r1,k1,alpha12,c1,mu,a,r2,alpha21,v,b,s,m,rho,gamma,delta,c2,lamb,c]

# Setup Time
tB = 0  #Beginning Time
tF = 2500 # Final Time
nTime = 10000 # time points

# Create a vector of time points from beginning time to end time that includes nTime number of points
tspan = np.linspace(tB,tF,nTime)

ode_out = odeint(model, y0, tspan ,args=(p,), hmax=1/100)

N1 = ode_out[:,0]
N2 = ode_out[:,1]
I = ode_out[:,2]
Q = ode_out[:,3]




# Plot Results
fig, (ax) = plt.subplots(1, 1)
ax.plot(tspan, I, 'tab:gray', label="I")  # Note you need the labels for the plot legend
ax.plot(tspan, N1, 'r', label="N1")
ax.plot(tspan, N2, 'b', label="N2")
ax.set_yscale('log')
ax.set_ylabel('Infected Cells (i)') # Vertical axis label
ax.set_xlabel('Time (days)') # Horizontal Axis label
ax.legend(loc='best') # Show legend and put it in the "best" location
ax.set_ylim(1e0,1e12+1)


# Zoom In
fig, (ax1,ax2) = plt.subplots(2, 1)
ax1.plot(tspan, I, 'tab:gray', label="I")  # Note you need the labels for the plot legend
ax1.plot(tspan, N1, 'r', label="N1")
ax1.plot(tspan, N2, 'b', label="N2")
ax1.set_yscale('log')
ax1.set_ylabel('Infected Cells (i)') # Vertical axis label
ax1.set_xlabel('Time (days)') # Horizontal Axis label
ax1.set_xlim(0, 140)
ax1.set_ylim(1e6,1e11)

# Q Function
ax2.plot(tspan,Q,label="Q")
ax2.set_ylabel('Q')
ax2.set_xlabel('Time (days)')
ax2.set_xlim(0,200)

plt.show() # Show Plot
