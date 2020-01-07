import numpy as np
from scipy.optimize import curve_fit


def kww(t, talpha, A, beta):
    return A*np.exp(-(t/talpha)**beta)


def kww_fit(xdata, ydata):
    param_bounds = ([np.finfo(float).eps, -np.inf, -np.inf],
                    [np.inf, np.inf, np.inf])
    popt, _ = curve_fit(kww, xdata, ydata, bounds=param_bounds)

    return popt
