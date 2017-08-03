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
 affiliations:
 - name: Institut Pasteur
   index: 1
 date: 2 August 2017
bibliography: paper.bib
---

# Summary

**Spectrum** is a Python library that includes tools to estimate Power Spectral Densities. Methods
available are based on Fourier transform, parametric methods or eigenvalues analysis. 

- The Fourier methods are based upon correlogram, periodogram and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are available (DPSS, Taylor, ...) [welch:1967],[marple:1987].
- The parametric methods are based on Yule-Walker, BURG, MA and ARMA, covariance and modified covariance methods [marple:1987],[percival:1993].
- Non-parametric methods based on eigen analysis (e.g., MUSIC) and minimum variance analysis are also implemented [marple:1987].
- Multitapering is also available [percival:1993]
- Classical tools useful to spectral analysis and more generally signal processing such as window tapering [harris:1978] or transfer function are also available within the library.

-![https://doi.org/10.6084/m9.figshare.5270866.v1](psd_all.png)

# References
