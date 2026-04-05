# ZPeak_CMS
Reconstruction and fitting of the dilepton invariant‑mass spectrum around the Z boson peak using CMS Open Data.
The project extracts muon/electron pairs, builds the invariant‑mass distribution, and fits the Z peak using Bayesian methods.
## Overview
This repository performs:
- Extraction of dilepton events from CMS AOD files
- Invariant‑mass reconstruction
- Modeling of the Z peak (Breit–Wigner + detector resolution)
- Parameter estimation using Dynesty and MCMC
- Generation of comparison plots and posterior distributions
## Project Files
|File|Description|
|----|-----------|
|file_extraction.py|Reads CMS Open Data and builds the invariant‑mass histogram|
|reconstruction.py|Reconstructs dilepton invariant mass from 4‑vectors|
|model.py|Defines the Z‑peak model and likelihood|
|main.ipynb|Main notebook for nested sampling fit + corner plot + comparison plot|
|mcmc_fitting.py|MCMC‑based parameter estimation|
|dinesty_fitting.py||Dynesty-based parameter esimation|
|Project Structure.pdf|Project instructions|
## Requirements
Install dependencies:
'pip install uproot numpy awkward scipy matplotlib uproot dynesty emcee corner'.

The CMS AOD data files are not stored in this repository because they are too large.
The required dataset must be downloaded from the CERN Open Data Portal and the file path must be manually updated inside 'file_extraction.py'.
## How to Run
1. Download the file from the CERN open data record. Since the files are quite big, they are not included in this repo. For this implementation we have decided to go with Mu_PAT_data_500files_01.root. In case the user wants to use a different file, download a different file and change the variable file_name in file_extraction.py.
2. Run the Jupyter Notebook main.ipynb. It does baysian statistical inference with both the emcee (Monte Carlo Markov Chain) and the dynesty (Nested Sampling) packages. For every library, it produces a corner plot relating the parameter correlation, a trace plot and a plot comparing the fitted model to the data.

Outputs: 
- Corner Plot with parameter correlation for both emcee and dynesty packages
- Trace Plot for both emcee and dynesty packages 
- Plot comparing fitted model with data for both emcee and dynesty packages

## Results
The repository includes:
- Posterior corner and trace plots for both inference methods
- Model vs. data comparison
- Extracted Z‑peak mass and width




