# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 20:06:33 2020

@author: harold romero
"""
from scipy.integrate import odeint
import matplotlib.pyplot as plt
"""modelo SEIRD   """
def derivSt( y, t, beta, gamma, sigma, mu):
    S, E, I, R, D = y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I - mu*I
    dRdt = gamma * I
    dDdt = mu*I
    return  dSdt, dEdt, dIdt, dRdt, dDdt

def SEIRD(S0, I0, beta, gamma, sigma,mu, days):
    y0 = S0,0, I0, 0,0
    t = list(range(0, days))
    result = odeint(derivSt, y0, t, args=(beta, gamma,sigma, mu))
    S, E, I, R, D = result.T
    return S, E, I, R, D
    
N      = 6700000
i_0    = 1
e_0    = 0
r_0    = 0
d_0    = 0
s_0    = N - i_0 - r_0
S, E, I, R,D = SEIRD(s_0, I0=i_0, beta=0.000001, gamma=1/2.9, sigma=1/5.2, mu=0.02, days=100)


""" modelo SIR """
def derivSti( y, t, betai, gammai):
    Si, Ii , Ri= y
    dSidt = -betai * Si * Ii
    dIidt= betai * Si * Ii - gammai * Ii
    dRidt = gammai * Ii
    return  dSidt, dIidt,dRidt

def SIR(S0, I0, betai, gammai, days):
    # Initial conditions vector
    y0 = S0, I0, 0
    # Integrating the SIR equations over the time grid, t
    t = list(range(0, days))
    # Getting results
    result = odeint(derivSti, y0, t, args=(betai, gammai))
    Si, Ii, Ri = result.T
    return Si, Ii, Ri

Si, Ii, Ri = SIR(s_0, I0=i_0, betai=0.000001, gammai=1/5.2, days=65)

"""modelo SEIR """
def derivSt_SEIR( y, t, betaii, gammaii,sigmaii):
    Sii, Eii, Iii, Rii = y
    dSiidt = -betaii * Sii * Iii
    dEiidt= betaii*Sii*Iii - sigmaii*Eii
    dIiidt= sigmaii*Eii - gammaii * Iii 
    dRiidt = gammaii * Iii
    return  dSiidt,dEiidt, dIiidt, dRiidt

def SEIR(S0, I0, betaii, gammaii, sigmaii, days):
    y0 = S0, 0, I0, 0
    t = list(range(0, days))
    result = odeint(derivSt_SEIR, y0, t, args=(betaii, gammaii, sigmaii))
    Sii, Eii, Iii, Rii = result.T
    return  Sii, Eii, Iii, Rii
    
Sii, Eii, Iii, Rii = SEIR(s_0, I0=i_0, betaii=0.000001, gammaii=1/2.9, sigmaii=1/5.2, days=65)
    
""" graficando"""
plt.plot(S)
plt.plot(E)
plt.plot(I)
plt.plot(R)
plt.plot(D)
plt.plot(Si)
plt.plot(Ii)
plt.plot(Ri)
plt.plot(Sii)
plt.plot(Eii)
plt.plot(Iii)
plt.plot(Rii)
plt.legend (['s','e','i','r','Si', 'Ii', 'Ri','Sii', 'Iii', 'Rii'])
plt.grid()
plt.show ()
