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

- The Fourier methods are based upon correlogram, periodogram and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are available (DPSS, Taylor, ...).
- The parametric methods are based on Yule-Walker, BURG, MA and ARMA, covariance and modified covariance methods.
- Non-parametric methods based on eigen analysis (e.g., MUSIC) and minimum variance analysis are also implemented.
- Multitapering is also available

-![Fidgit deposited in figshare.](figshare_article.png)

# References
  
