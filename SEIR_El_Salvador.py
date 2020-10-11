import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate, optimize

np.set_printoptions(suppress=True)

datos = pd.read_csv('El Salvadr_Activos.csv')
df=pd.DataFrame(datos)
nmp=df.to_numpy()
activos=(nmp[3:4])
recuperados= (nmp[1:2])
j , dias = activos.shape 
ejex = list(range(0, dias))
activos = activos[0]

recuperados=recuperados[0]

#activos = np.array(activos, dtype=float)
#ejex = np.array(ejex, dtype=float)
#plt.plot(activos)

      
def derivSt_SIR( y, t, beta, gamma):
    S,  I, R = y
    dSdt = -beta * S * I
    dIdt=  beta * S * I - gamma * I 
    dRdt = gamma * I
    return  dSdt, dIdt, dRdt

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

def sir_model(y, x, beta, gamma):
    S = -beta * y[0] * y[1] / N
    R = gamma * y[1]
    I = -(S + R)
    return S, I, R

def fit_odeint(x, beta, gamma):
    return integrate.odeint(sir_model, (S0, I0, R0), x, args=(beta, gamma))[:,1]


ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
xdata = list(range(0, 14))

#ydata= activos
#ydata1=recuperados
#xdata= ejex

#ydata = np.array(ydata, dtype=float)
#ydata1= np.array(ydata1, dtype=float)
#xdata = np.array(xdata, dtype=float)

N = 763
I0 = ydata[0]
S0 = N - I0
R0 = 0

popt, pcov = optimize.curve_fit(fit_odeint, xdata, ydata)
#pop , pcov = optimize.curve_fit(fit_odeint, xdata, ydata1)
fitted = fit_odeint(xdata, *popt)
#fitted1= fit_odeint(xdata, *pop)
beta1, gamma1 = popt

SIR(S0, I0,  beta1/N, gamma1, 300)


