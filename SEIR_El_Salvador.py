#UNIVERSIDAD DE EL SALVADOR
#FACULTAD DE INGENIERIA Y ARQUITECTURA
#ESCUELA DE INGENIERIA EECTRICA
#PROYECTO DE INGENIERIA ELECTRICA 1
#DOCENTE ASESOR: PHD. CARLOS MARTINEZ
#ESTUDIANTES: BR. FIDEL ERNESTO CASTRO CONTRERAS CC14009
#             BR. HAROLD ERNESTO ROMERO PALACIOS RP12010
             
# EL PRESENTE CODIGO TIENE LA FUNCION DE CARGAR LOS DATOS DE CASOS CONFIRMADOS
# DE COVID-19 EN EL SALVADOR, PARAMETRIZARLOS EN MODELO SIR UTILIZANDO EL 
# METODO DE MINIMOS CUADRADOS Y OBTENER PROYECCIONES PARA LOS MODELOS SIR,
# SEIR Y SEIRD PARA 300 DIAS DESPUES DEL PRIMER CONTAGIO EN EL SALVADOR

### SE IMPORTAN LAS LIBRERIAS CON LAS QUE SE TRABAJARAN

import pandas as pd 
# SE UTILIZARA PARA CARGAR LOS DATOS .CSV
import numpy as np
import matplotlib.pyplot as plt 
#SE UTILIZA PARA GRAFICAR LOS RESULTADOS
from scipy.integrate import odeint 
#SE UTILIZA PARA RESOLVER LAS ECUACIONES DIFERENCIALES
from scipy import integrate, optimize
#SE UTILIZA PARA AJUSTAR LOS DATOS Y OBTENER LOS PARAMETROS DEL MODELO SIR

#DEFINICION DE LOS MODULOS A UTLIZAR

def derivSt_SIR( y, t, beta, gamma):
    #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
    #SIR, RETORNA LAS DERIVADAS
    S,  I, R = y
    dSdt = -beta * S * I
    dIdt=  beta * S * I - gamma * I 
    dRdt = gamma * I
    return  dSdt, dIdt, dRdt

def derivSt_SEIR( y, t, beta, gamma,sigma=1/5.2):
    S, E, I, R = y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I 
    dRdt = gamma * I
    return  dSdt,dEdt, dIdt, dRdt

def derivSt_SEIRD( y, t, beta, gamma, sigma=1/5.2, mu=0.02):
    S,  E, I, R, D= y
    dSdt = -beta * S * I
    dEdt= beta*S*I - sigma*E
    dIdt= sigma*E - gamma * I - mu*I
    dRdt = gamma * I
    dDdt = mu*I
    return  dSdt, dEdt, dIdt, dRdt, dDdt

def SIR(S0, I0, beta, gamma, days):
    #MODULO UTILIZADO PARA RESOLVER LAS ECUACIONES DIFERENCIALES DEL MODELO
    #SIR, RECIBE CANTIDAD DE SUCEPTIBLES INICIALES (POBLACION DEL PAIS), 
    #INFECTADOS INICIALES, PARAMETROS BETA Y GAMMA DEL MODELO SIR Y LA 
    #CANTIDAD DE DIAS A LOS QUE SE QUIERE HACER LA PROYECCION
    y0 = S0, I0, 0
    t = list(range(0, days))
    result = odeint(derivSt_SIR, y0, t, args=(beta, gamma))
    S, I, R = result.T
    plt.plot(S)
    plt.plot(I)
    plt.plot(R)
    #plt.plot(ydata, 'o')
    plt.grid()
    plt.legend(['s',"i",'r'])
    #plt.legend(['datos ajustados', 'datos sin procesar'])
    # plt.title("Modelo SIR")
    plt.title("PROYECCION COVID-19 EL SALVADOR PARA 300 DIAS DESPUES DEL PRIMER CONTAGIO MODELO SIR")
    plt.show()
    return 0

def SEIR(S0, I0, beta, gamma,  days):
    y0 = S0, 0, I0, 0
    t = list(range(0, days))
    result = odeint(derivSt_SEIR, y0, t, args=(beta, gamma))
    S, E, I, R = result.T
    plt.plot(S)
    plt.plot(E)
    plt.plot(I)
    plt.plot(R)
    plt.grid()
    plt.legend(['s',"e","i",'r'])
    plt.title("PROYECCION COVID-19 EL SALVADOR PARA 300 DIAS DESPUES DEL PRIMER CONTAGIO MODELO SEIR")
    plt.show()
    return 0
    
def SEIRD(S0, I0, beta, gamma, days):
    y0 = S0,0, I0, 0,0
    t = list(range(0, days))
    result = odeint(derivSt_SEIRD, y0, t, args=(beta, gamma))
    S, E, I, R, D = result.T
    plt.plot(S)
    plt.plot(E)
    plt.plot(I)
    plt.plot(R)
    plt.plot(D)
    plt.grid()
    plt.legend(['s',"e","i",'r',"d"])
    plt.title("PROYECCION COVID-19 EL SALVADOR PARA 300 DIAS DESPUES DEL PRIMER CONTAGIO MODELO SEIRD")
    plt.show()
    return 0

def sir_model(y, x, beta, gamma):
    #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
    #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    #y[0] ES EL VECTOR CON SUCEPTIBLES Y y[1] ES EL VECTOR CON INFECTADOS
    S = -beta * y[0] * y[1] / N
    R = gamma * y[1]
    I = -(S + R)
    return S, I, R

def seir_model(y, x, beta, gamma):
    #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
    #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    #y[0] ES EL VECTOR CON SUCEPTIBLES Y y[1] ES EL VECTOR CON EXPUESTOS Y y[2]
    #VECTOR DE INFECTADOS 
    S = -beta * y[0] * y[2] / N
    E = (S + 1/5.2*y[1])
    R = gamma * y[2]
    I = 1/5.2*y[1]-R
    return S,E, I, R

def seird_model(y, x, beta, gamma, sigma =1/5.2, mu=0.02):
    #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
    #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    #y[0] ES EL VECTOR CON SUCEPTIBLES Y y[1] ES EL VECTOR CON EXPUESTOS Y y[2]
    #VECTOR DE INFECTADOS 
    S = -beta * y[0] * y[2] / N
    E = (S + sigma*y[1])
    R = gamma * y[2]
    I = sigma*y[1]-R - mu*y[2]
    D = mu*y[2]
    return S,E, I, R, D


def fit_odeint_SIR(x, beta, gamma):
     #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
     #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    return integrate.odeint(sir_model, (S0, I0, R0), x, args=(beta, gamma))[:,1]

def fit_odeint_SEIR(x, beta, gamma):
     #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
     #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    return integrate.odeint(seir_model, (S0, E0, I0, R0), x, args=(beta, gamma))[:,1]

def fit_odeint_SEIRD(x, beta, gamma, sigma, mu):
     #MODULO UTILIZADO PARA DESCRIBIR LAS ECUACIONES DIFERENCIALES DEL MODELO
     #SIR QUE SE UTILIZARAN EN PARA EL AJUSTE DE DATOS POR MINIMOS CUADRADOS
    return integrate.odeint(seird_model, (S0, E0, I0, R0, D0), x, args=(beta, gamma))[:,1]

def ESIR():
    popt, pcov = optimize.curve_fit(fit_odeint_SIR, xdata, ydata)
    #RETORNA LOS PARAMETROS BETA Y GAMMA OBTENIDOS POR EL METODO DE MINIMOS 
    #CUADRADOS
    fitted = fit_odeint_SIR(xdata, *popt)
    #SE AJUSTA LA CURVA CON LOS PARAMETROS CONTENIDOS EN "popt"
    beta, gamma = popt
    #SE GUARDAN LOS PARAMETROS BETA Y GAMMA
    print(('Beta = ', beta, 'Gamma = ', gamma))
    #SE IMPRIME EN CONSOLA BETA Y GAMMA
    SIR(S0, I0,  beta/N, gamma, 300)
    return 0

def ESEIR():
    popt, pcov = optimize.curve_fit(fit_odeint_SEIR, xdata, ydata)
    fitted = fit_odeint_SEIR(xdata, *popt)
    beta, gamma = popt
    print(('Beta = ', beta, 'Gamma = ', gamma))
    SEIR(S0, I0,  beta/N, gamma,300)
    return 0

def ESEIRD():
    popt, pcov = optimize.curve_fit(fit_odeint_SEIRD, xdata, ydata)
    fitted = fit_odeint_SEIRD(xdata, *popt)
    beta, gamma,sigma,mu = popt
    print(('Beta = ', beta, 'Gamma = ', gamma))
    SEIRD(S0, I0,  beta/N, gamma, 300)
    return 0

np.set_printoptions(suppress=True)
#SE CANCELA EL USO DE NOTACION CIENTIFICA EN LA VISUALIZACION DE LOS RESULTADOS

datos = pd.read_csv('El Salvadr_Activos.csv')
#SE CARGAN LOS DATOS DE EL SALVADOR
df=pd.DataFrame(datos)
#SE CONVIENTERN LOS DATOS TRAIDOS CON PANDAS A UN DATAFRAME
nmp=df.to_numpy()
#SE CONVIERTE EL DATAFRAME A UNA MATRIZ NUMPY
activos=(nmp[3:4])
#SE CONSERVAN UNICAMENTE LOS CASOS ACTIVOS
recuperados= (nmp[1:2])
#SE CONSERVAN UNICAMENTE LOS CASOS RECUPERADOS
j , dias = activos.shape 
#SE CGUARDAN LA CANTIDAD DE DIAS TRASCURRIDOS DESDE EL PRIMER CONTAGIO
ejex = list(range(0, dias))
#SE CREA EL EJE X CON EL QUE SE GRAFICARA
activos = activos[0]
#PARA SIMPLICAR EL PROCESO SE DA LA INSTRUCCION DE CONSERVAR EL VECTOR 0 DE LA
#MATRIZ DE ACTIVOS
recuperados=recuperados[0]
#PARA SIMPLICAR EL PROCESO SE DA LA INSTRUCCION DE CONSERVAR EL VECTOR 0 DE LA
#MATRIZ DE RECUPERADOS



ydata= activos
#SE PASA EL VECTOR CON CASOS ACTIVOS AL VECTOR "ydata"
xdata= ejex
#SE PASA EL VECTOR CON LA CANTIDAD DEL DIAS AL VECTOR "ejex"

ydata = np.array(ydata, dtype=float)
xdata = np.array(xdata, dtype=float)
#SE FORZA LA CONVERSION DE LOS DATOS CONTENIDOS EN "ydata" y "xdata" al tipo
#FLOTANTE

N = 67000000
#POBLACION DE EL SALVADOR
I0 = 1
#NUMERO DE INFECTADOS INICIALES
E0=0
#NUMERO DE EXPUESTOS INICIALES
S0 = N - I0
#NUMERO DE SUCEPTIBLES
R0 = 0
#CANTIDAD DE RECUPERADOS INICIALES
D0 = 0
#CANTIDAD DE FALLECIDOS INICIALES

ESEIRD()
#SE EJECUTA EL MODELO SEIR




#plt.title("Influenza in a boarding school (British Medical Journal, 4 March 1978")
#ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
#xdata = list(range(0, 14))






