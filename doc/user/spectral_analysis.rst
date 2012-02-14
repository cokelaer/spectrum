
Quick overview of spectral analysis methods
###############################################

This section gives you a quick overview of the spectral analysis methods and classes that are available in **spectrum**. You will find the different classes associated to each PSD estimates. A functional approach is also possible but is not described here. See the reference guide for more details.


Non-parametric classes
=========================

The Fourier-based methods provides :class:`~spectrum.periodogram.Periodogram`, :class:`~spectrum.correlog.pcorrelogram`, Welch estimate (not implemented see pylab.psd instead) and multitapering (not yet implemented).

In addition to the Fourier-based methods, there are 3 types of non-parametric methods:

#. The Minimum of variance MV (Capon) is implemented in the class :class:`~spectrum.minvar.pminvar`.
#. Two eigenvalues decomposition (MUSIC, eigenvalue) can be found in :class:`~spectrum.eigenfreq.pev` and :class:`~spectrum.eigenfre.pmusic`.
#. Maximum entropy (MEM) (not yet implemented)

Autoregressive spectral estimation
========================================

There are essentially 3 methods to estimate the autoregressive (AR) parameters.
The first one uses the autocorrelation sequence such as in the so-called 
**Yule-Walker** method (see :class:`~spectrum.yulewalker.pyule`). A second method uses the reflection coefficient method such as in the **Burg** algorithm (see :class:`~spectrum.burg.pburg`). These methods minimise the forward prediction error (and backward) using Levinson recursion. Finally, a third important category of AR parameter method is based on the least squares linear prediction, which can be further decomposed into 2 categories. One that separate the minimization of the forward and backward linear prediction squared errors such as the 
**autocorrelation** or **covariance** methods (see :class:`~spectrum.covar.pcovar`). Another one that performs a 
combined minimization of the forward and backward prediction squared errors 
(**modified covariance**) (see :class:`~spectrum.modcovar.pmodcovar`). 

Spectrum also provides :class:`~spectrum.arma.parma`, :class:`~spectrum.arma.pma` classes.

