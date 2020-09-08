
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize

ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
xdata = np.linspace(0,14, dtype="float")



ydata = np.array(ydata, dtype=float)
xdata = np.array(xdata, dtype=float)

def sir_model(y, x, beta, gamma):
    S = -beta * y[0] * y[1] / N
    R = gamma * y[1]
    I = -(S + R)
    return S, I, R

def fit_odeint(x, beta, gamma):
    return integrate.odeint(sir_model, (S0, I0, R0), x, args=(beta, gamma))[:,1]

N = sum(ydata)
I0 = ydata[0]
S0 = N - I0
R0 = 0.0

popt, pcov = optimize.curve_fit(fit_odeint, xdata, ydata)
fitted = fit_odeint(xdata, *popt)

plt.plot(xdata, ydata, 'o')
plt.plot(xdata, fitted)
plt.show()
