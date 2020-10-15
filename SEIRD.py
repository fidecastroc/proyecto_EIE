from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.special import logit, expit
import matplotlib.pyplot as plt
import numpy as np

def derivSt_SIR( y, t, beta, gamma):
    S,  I, R = y
    dSdt = -beta * S * I
    dIdt=  beta * S * I - gamma * I 
    dRdt = gamma * I
    return  dSdt, dIdt, dRdt

def derivSt_SEIR( y, t, beta, gamma,sigma):
    S, E, I, R = y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I 
    dRdt = gamma * I
    return  dSdt,dEdt, dIdt, dRdt

def derivSt_SEIRD( y, t, beta, gamma, sigma, mu):
    S,  E, I, R, D= y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I - mu*I
    dRdt = gamma * I
    dDdt = mu*I
    return  dSdt, dEdt, dIdt, dRdt, dDdt

def SIR(S0, I0, beta, gamma, days):
    y0 = S0, I0, 0
    t = list(range(0, days))
    result = odeint(derivSt_SIR, y0, t, args=(beta, gamma))
    S, I, R = result.T
    plt.plot(S)
    plt.plot(I)
    plt.plot(R)
    plt.grid()
    plt.legend(['s',"i",'r'])
    plt.title("Modelo SIR")
    plt.show()
    return 0

def SEIR(S0, I0, beta, gamma, sigma, days):
    y0 = S0, 0, I0, 0
    t = list(range(0, days))
    result = odeint(derivSt_SEIR, y0, t, args=(beta, gamma, sigma))
    S, E, I, R = result.T
    plt.plot(S)
    plt.plot(E)
    plt.plot(I)
    plt.plot(R)
    plt.grid()
    plt.legend(['s',"e","i",'r'])
    plt.title("Modelo SEIR")
    plt.show()
    return 0
    
def SEIRD(S0, I0, beta, gamma, sigma,mu, days):
    y0 = S0,0, I0, 0,0
    t = list(range(0, days))
    result = odeint(derivSt_SEIRD, y0, t, args=(beta, gamma, sigma, mu))
    S, E, I, R, D = result.T
    plt.plot(S)
    plt.plot(E)
    plt.plot(I)
    plt.plot(R)
    plt.plot(D)
    plt.grid()
    plt.legend(['s',"e","i",'r',"d"])
    plt.title("Modelo SEIRD")
    plt.show()
    return 0
    
def SIR_SSQ(coeff, N, y0, Iob):
    beta, gamma = coeff
    t = list (range(0,len(Iob)))
    result = odeint(derivSt_SIR, y0, t, args=(beta, gamma))
    S,I,R = result.T
    f = np.sum((I-Iob)**2)
    return f

N      = 6700000
i_0    = 1
e_0    = 0
r_0    = 0
d_0    = 0
s_0    = N - i_0 - r_0


S(s_0, I0=i_0, beta=3.6029094743522356/N, gamma=1/3.533526760861618,sigma=1/5.2,mu=0.02,days=300)


