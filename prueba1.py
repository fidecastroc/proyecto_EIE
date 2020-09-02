from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import logit, expit
import matplotlib.pyplot as plt
import numpy as np

#SIR

def derivSt( y, t, beta, gamma, sigma, mu):
    S, E, I, R, D = y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I - mu*I
    dRdt = gamma * I
    dDdt = mu*I
    return  dSdt, dEdt, dIdt, dRdt, dDdt

def SEIRD(S0, I0, beta, gamma, sigma,mu, days):
    # Initial conditions vector
    y0 = S0,0, I0, 0,0
    # Integrating the SIR equations over the time grid, t
    t = list(range(0, days))
    # Getting results
    result = odeint(derivSt, y0, t, args=(beta, gamma,sigma, mu))
    S, E, I, R, D = result.T
    return S, E, I, R, D
    
     
N = 6700000
i_0    = 1
e_0    = 0
r_0    = 0
d_0    = 0
s_0    = N - i_0 - r_0
S, E, I, R,D = SEIRD(s_0, I0=i_0, beta=0.000001, gamma=1/2.9, sigma=1/5.2, mu=0.02, days=100)

plt.plot(S)
plt.plot(E)
plt.plot(I)
plt.plot(R)
plt.plot(D)
plt.grid()
plt.legend(["s","e","i","r","d"])
plt.show()
