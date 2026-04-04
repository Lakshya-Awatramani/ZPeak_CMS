## Importing Libraries
import numpy as np 
import dynesty
from reconstruction import m_inv_mass2
from model import log_likelihood
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt

## Importing data
# Invaraint mass data
data = np.sqrt(m_inv_mass2)

bin_edges = np.linspace(70, 100, 500) 
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Number of counts in a bin
observed_counts, _ = np.histogram(data, bins=bin_edges)

## Prepare for fitting with dynesty
ndim = 5

# Define 5D Correlation Matrix 
C = np.identity(ndim)
Cinv = np.linalg.inv(C)

#Define a uniform prior transform
def prior_transform(u):
    '''
    Maps a unit cube to the parameter space of the model.
     - m0: mass of the Z boson, expected around 91 GeV,
     - gamma: width of the Z boson, expected around 2.5 GeV,
     - sigma: width of the Voigt profile, expected around 1 GeV,
     - Nsig: total signal events,
     - Nbg: background level per bin.
    The parameters are scaled to reasonable ranges based on the observed data.

    Parameters:
    u : array-like
        A 5D array of values in the unit cube [0,1]^5.
    Returns:
    array: A 5D array of model parameters [m0, gamma, sigma, Nsig, Nbg].
    '''
    u = np.asarray(u)
    # ensure observed_counts works if it's an awkward array
    obs = np.asarray(observed_counts, dtype=float)
    total_counts = float(np.sum(obs))
    # background per-bin upper bound: use either max observed per-bin or average
    max_bg_per_bin = float(max(obs.max(), total_counts / len(bin_centers)))
    eps = 1e-8
    m0 = 80.0 + 20.0 * u[0]        # maps [0,1] -> (80,100)
    gamma = 5.0 * u[1]             # maps [0,1] -> (0,5)
    sigma = 5.0 * u[2]             # maps [0,1] -> (0,5)
    Nsig = eps + u[3] * (total_counts - eps)      # total signal events
    Nbg = eps + u[4] * (max_bg_per_bin - eps)     # background level per bin
    return np.array([m0, gamma, sigma, Nsig, Nbg])

# Use the log_likelihood function with the prior_transform
sampler = dynesty.NestedSampler(log_likelihood, prior_transform, ndim,
                                logl_args=(bin_centers, observed_counts))
sampler.run_nested()
results = sampler.results