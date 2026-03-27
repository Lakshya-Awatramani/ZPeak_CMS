import scipy as sc
import numpy as np

def model(x, A, sigma, gamma, bg_val):
    """
    Voigt profile (Breit-Wigner ⊗ Gaussian) with amplitude,
    plus a linear background.
      A         - peak amplitude
      sigma     - Gaussian std dev (resolution)
      gamma     - Lorentzian HWHM (natural/decay width)
      bg_val    - background intercept
    """
    signal = A * sc.special.voigt_profile(x, sigma, gamma)
    background = bg_val
    return signal + background

def log_prior(C):
    """Uniform prior with physical boundary enforcement."""
    A, sigma, gamma, bg_val, bg_slope = C
    if A > 0 and sigma > 0 and gamma > 0:
        return 0.0
    return -np.inf

def log_likelihood(data, sigma_data, C):
    """Gaussian log-likelihood with per-point uncertainties."""
    residuals = (data - model(x, *C)) / sigma_data
    return -0.5 * np.sum(residuals**2)

def log_posterior(C, data, sigma_data):
    lp = log_prior(C)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(data, sigma_data, C)