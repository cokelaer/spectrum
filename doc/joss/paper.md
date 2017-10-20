---
title: 'Spectrum': Spectral Analysis in Python
tags:
  - spectral analysis
  - periodogram
  - yule-walker
  - multi-tapering
  - burg
  - ARMA
  - eigen-values
  - tapering windows
authors:
 - name: Thomas Cokelaer
   orcid: 0000-0001-6286-1138
   affiliation: 1
 - name: Juergen Hasch
   orcid: 0000-0002-9457-1220
   affiliation: 2
 affiliations:
 - name: Institut Pasteur - Bioinformatics and Biostatistics Hub - C3BI, USR 3756 IP
CNRS - Paris, France
   index: 1
 - name:  Robert Bosch GmbH: Renningen, Baden-Württemberg, Germany
   index: 2
date: 2 August 2017
bibliography: paper.bib
---

# Summary

**Spectrum** is a Python library that includes tools to estimate Power Spectral Densities. Although the use of 
power spectrum of a signal is fundamental in electrical engineering (e.g. radio communications, radar), it has
a wide range of applications from cosmology (e.g., detection of gravitational waves in 2016), to music 
(pattern detection) or biology (mass spectroscopy).

Methods available are based on Fourier transform, parametric methods or eigenvalues analysis. Although standard methods such as 
periodogram are available, less common methods (e.g. multitapering) are also implemented:

- The Fourier methods are based upon correlogram, periodogram and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are available (DPSS, Taylor, ...)[harris:1978],[welch:1967],[marple:1987].
- The parametric methods are based on Yule-Walker, BURG, MA and ARMA, covariance and modified covariance methods [marple:1987],[percival:1993].
- Non-parametric methods based on eigen analysis (e.g., MUSIC) and minimum variance analysis are also implemented [marple:1987].
- Multitapering method is also available [percival:1993]
- Classical tools useful to spectral analysis and more generally signal processing such as window tapering [harris:1978] or transfer function are also available within the library.

The following image shows the different methods of spectral estimation that are available in **Spectrum**.

-![https://doi.org/10.6084/m9.figshare.5270866.v1](psd_all.png)

**Spectrum** relies on Matplotlib [matplotlib:2007] for the plotting. We also
use Numpy [numpy:2011] for fast array manipulation and Scipy [scipy:2001] for
linear algebra.

# References
