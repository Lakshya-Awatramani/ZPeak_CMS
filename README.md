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
|Defines the Z‑peak model and likelihood|Nested sampling fit + corner plot + comparison plot|
|mcmc_fitting.py|MCMC‑based parameter estimation|
|bestfit_parameters.txt|Saved best‑fit values from Dynesty|
|Project Structure.pdf|Project instructions|
## Requirements
Install dependencies:
'pip install numpy scipy matplotlib uproot dynesty emcee'.

The CMS AOD data files are not stored in this repository because they are too large.
The required dataset must be downloaded from the CERN Open Data Portal and the file path must be manually updated inside 'file_extraction.py'.
## How to Run
1. Extract events:
'python file_extraction.py'
2. Reconstruct invariant mass:
'python reconstruction.py'
3. Fit with Dynesty: 
'python dynesty_fitting.py'

Outputs: 
- CornerPlotDynestyFitting.png
- Dynesty comparison plot.png
- bestfit_parameters.txt
4. (Optional) MCMC fit:
  'python mcmc_fitting.py'
## Results
The repository includes:
- Posterior corner plots
- Model vs. data comparison
- Extracted Z‑peak mass and width




