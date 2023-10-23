import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt 


# Define your ODE Model 
# y are state variables (dependent on time and parameters)
# t is for time
# k is for parameters (rates, proportions, etc.)
def model(y,t,k):  
    # State Variables
    S = y[0]
    I = y[1]
    C = y[2]
    W = y [3]



    # Parameters (changing variables)
    beta = k[0]
    Kone = k[1]
    K1 = k[2]
    K2 = k[3]
    A = k[4]
    A0 = k[5]
    KK = k[6]
    K0 = k[7]
    K1 = k[8]
    B1 = k[9]
    B0 = k[10]
    b = k[11]
    BB = k[12]
    beta0 = k[13]
    beta1 = k[14]
    
    beta,Kone, K1, K2, A, A0, K, K0, K1, B1, B0, b, BB, beta0, beta1


    # Differential Equations
    dsdt = A-A0*S-(beta*S*C)
    didt = beta*S*C - beta0*I - beta1*C*I
    dcdt = KK-K0*C-K1*C*W
    dwdt = BB+b*C-B0*W-B1*W*C
    
    
    # Compile and return    
    dydt = [dsdt, didt, dcdt, dwdt]
    return dydt

# Initial Conditions
S0 = 149
I0 = 1.26
C0 = 1.95
W0 = 252.26


y0 = [S0, I0, C0, W0]



# Parameter Values
beta = 0.00005
Kone = 10
K2 = 11
K3 = 12
A = 1.5
A0 = 0.01
K = 12
K0 = 5 
K1 = 0.005
B1 = 0.001
B0 = 0.05
b = 0.01
BB = 15
beta0 = 0.003
beta1 = 0.005

p = [beta,Kone, K1, K2, A, A0, K, K0, K1, B1, B0, b, BB, beta0, beta1]

# Setup Time
tB = 0  #Beginning Time
tF = 100 # Final Time
nTime = 100 # time points

# Create a vector of time points from beginning time to end time that includes nTime number of points
tspan = np.linspace(tB,tF,nTime)

ode_out = odeint(model, y0, tspan, args=(p,))

S = ode_out[:,0]
I = ode_out[:,1]
C = ode_out[:,2]
W = ode_out[:,3]



# Setup Time
tB = 0  #Beginning Time
tF = 100 # Final Time
nTime = 100 # time points



# Plot Results
plt.plot(tspan, C, 'r:', label="Infectious")  # Note you need the labels for the plot legend
plt.ylabel('Infected Cells (i)') # Vertical axis label
plt.xlabel('Time (days)') # Horizontal Axis label
plt.legend(loc='best') # Show legend and put it in the "best" location
plt.show() # Show PLot



