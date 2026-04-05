import numpy as np
from scipy.special import voigt_profile


def voigt(x, m0, gamma, sigma):
    return voigt_profile(x - m0, sigma, gamma)

def model(bin_centers, m0, gamma, sigma, Nsig, Nbg):
    """Predicted counts per bin: signal + flat background"""
    signal = voigt(bin_centers, m0, gamma, sigma)
    signal *= Nsig / np.sum(signal)  # scale to total signal events
    background = np.full_like(bin_centers, Nbg)
    return signal + background

def log_likelihood(theta, bin_centers, observed_counts):

    """Poisson log likelihood"""

    m0, gamma, sigma, Nsig, Nbg = theta
    mu = model(bin_centers, m0, gamma, sigma, Nsig, Nbg)
    
    return np.sum(observed_counts * np.log(mu) - mu)

def log_prior(theta):
    """Flat prior within bounds"""
    m0, gamma, sigma, Nsig, Nbg = theta
    if not (80 < m0 < 100 and 0 < gamma < 5 and 0 < sigma < 5 and Nsig > 0 and Nbg > 0):
        return -np.inf
    return 0  

def log_posterior(theta, bin_centers, observed_counts):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, bin_centers, observed_counts)