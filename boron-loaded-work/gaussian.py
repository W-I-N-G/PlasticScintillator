import numpy as np
from scipy.special import erfc
from scipy import optimize

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)


def gaussian_exp(x, amplitude, mean, stddev, a, b):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2) + a*np.exp(-b*x)


def gaussian_power(x, amplitude, mean, stddev, a, b):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2) + a*x**(-b)

def gauss_fit(x,y,params):
    popt, _ = optimize.curve_fit(gaussian, x, y, p0=params)
    u = popt[1]
    s = abs(popt[2])
    g = gaussian(x, *popt)
    return(g,u,s)

def gauss_fit_power(x,y,params):
    popt, _ = optimize.curve_fit(gaussian_power, x, y, p0=params)
    u = popt[1]
    s = abs(popt[2])
    g = gaussian_power(x, *popt)
    return(g,u,s)

def gauss_exp_fit(x,y,params):
    popt, _ = optimize.curve_fit(gaussian_exp, x, y,p0=params)
    u = popt[1]
    s = abs(popt[2])
    g = gaussian_exp(x, *popt)
    return(g,u,s)
