from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt

ydata = [3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4]
xdata = np.linspace(0,14)
yhat = savgol_filter(ydata, 13, 11)



plt.plot(yhat, color='red')
plt.show()
