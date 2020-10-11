import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate, optimize

<<<<<<< HEAD:untitled0.py
np.set_printoptions(suppress=True)

datos = pd.read_csv('El Salvadr_Activos.csv')
df=pd.DataFrame(datos)
nmp=df.to_numpy()
activos=(nmp[3:4])
j , dias = activos.shape 
ejex = list(range(0, dias))
activos = activos[0]
#activos = np.array(activos, dtype=float)
#ejex = np.array(ejex, dtype=float)
#plt.plot(activos)

      
=======
ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
xdata = list(range(0, 14))

ydata = np.array(ydata, dtype=float)
xdata = np.array(xdata, dtype=float)

N = 763
I0 = ydata[0]
S0 = N - I0
R0 = 0.0
"""" modelo sir"""
>>>>>>> 421a5fa106d4985326c8a698723b75d5777dca9a:ordenar.py
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
    #plt.plot(S)
    plt.plot(I)
    #plt.plot(R)
    plt.grid()
    plt.legend(['s',"i",'r'])
    plt.title("Modelo SIR")
    plt.show()
    return 0

"""sir para ajustar los nuevos beta y gamma """
def sir_model(y, x, beta, gamma):
    S = -beta * y[0] * y[1] / N
    R = gamma * y[1]
    I = -(S + R)
    return S, I, R

def fit_odeint(x, beta, gamma):
    return integrate.odeint(sir_model, (S0, I0, R0), x, args=(beta, gamma))[:,1]

<<<<<<< HEAD:untitled0.py

#ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
#xdata = list(range(0, 14))

ydata= activos
xdata= ejex

ydata = np.array(ydata, dtype=float)
xdata = np.array(xdata, dtype=float)

N = 6000000
I0 = ydata[0]
S0 = N - I0
R0 = 0

popt, pcov = optimize.curve_fit(fit_odeint, xdata, ydata)
fitted = fit_odeint(xdata, *popt)
beta1, gamma1 = popt
SIR(S0, I0,  beta1/N, gamma1, 400)

plt.plot(xdata, ydata, 'o')
plt.plot(xdata, fitted)
#plt.show()
=======
"""obteniendo los valores de beta y gamma que mejor se ajusten"""
popt, pcov = optimize.curve_fit(fit_odeint, xdata, ydata)
fitted = fit_odeint(xdata, *popt)
beta1, gamma1 = popt
print(beta1,gamma1)

""" nuevos parametros para el modelo"""
SIR(S0, I0,  beta1/N, gamma1,14)

""" graficando"""
plt.plot(xdata, ydata, 'ro')
plt.plot(xdata, fitted)
plt.show()
>>>>>>> 421a5fa106d4985326c8a698723b75d5777dca9a:ordenar.py
