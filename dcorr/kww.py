import numpy as np
from scipy.optimize import curve_fit


def kww(t, ta, alpha, beta):
    return alpha*np.exp(-(t/ta)**beta)


def kww_fit(xdata, ydata):
    popt, _ = curve_fit(kww, xdata, ydata)

    return popt
