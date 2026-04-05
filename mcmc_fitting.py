import numpy as np
import matplotlib.pyplot as plt
import emcee
from reconstruction import m_inv_mass2
from model import *

# Invaraint mass data
data = np.sqrt(m_inv_mass2)

bin_edges = np.linspace(70, 100, 500) 
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Number of counts in a bin
observed_counts, _ = np.histogram(data, bins=bin_edges)

ndim = 5
nwalkers = 50
nsteps = 3000

# Initial guesses near values seen in the histogram
starting_guesses = np.array([91, 2.0, 2.0, max(observed_counts)*10, np.median(observed_counts)]) + 1e-1 * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(bin_centers, observed_counts))
sampler.run_mcmc(starting_guesses, nsteps, progress=True)

chain = sampler.get_chain()
samples = sampler.get_chain(discard=500, flat=True)
best_fit = np.median(samples, axis=0)
