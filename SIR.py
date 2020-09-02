from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import logit, expit
import matplotlib.pyplot as plt
import numpy as np

#SIR

def derivSt( y, t, beta, gamma):
    S, I , R= y
    dSdt = -beta * S * I
    dIdt= beta * S * I
    dRdt = gamma * I
    return  dSdt, dIdt,dRdt

def SIR(S0, I0, beta, gamma, days):
    # Initial conditions vector
    y0 = S0, I0, 0
    # Integrating the SIR equations over the time grid, t
    t = list(range(0, days))
    # Getting results
    result = odeint(derivSt, y0, t, args=(beta, gamma))
    S, I, R = result.T
    return S, I, R
    
     
n = 6700000
i_0    = 1
r_0    = 0
s_0    = N - i_0 - r_0
S, I, R = SIR(s_0, I0=i_0, beta=0.00000002, gamma=1/5.2, days=165)

plt.plot(S)
plt.plot(I)
plt.plot(R)
plt.show()