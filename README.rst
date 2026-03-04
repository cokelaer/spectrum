SPECTRUM : Spectral Analysis in Python
==========================================

.. image:: https://badge.fury.io/py/spectrum.svg
    :target: https://pypi.python.org/pypi/spectrum

.. image:: https://github.com/cokelaer/spectrum/actions/workflows/main.yml/badge.svg?branch=main
    :target: https://github.com/cokelaer/spectrum/actions/workflows/main.yml

.. image:: https://coveralls.io/repos/cokelaer/spectrum/badge.png?branch=main
    :target: https://coveralls.io/r/cokelaer/spectrum?branch=main

.. image:: https://anaconda.org/conda-forge/spectrum/badges/license.svg
   :target: https://anaconda.org/conda-forge/spectrum

.. image:: https://anaconda.org/conda-forge/spectrum/badges/version.svg
    :target: https://anaconda.org/conda-forge/spectrum

.. image:: https://anaconda.org/conda-forge/spectrum/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/spectrum

.. image:: https://joss.theoj.org/papers/e4e34e78e4a670f2ca9a6a97ce9d3b8e/status.svg
   :target: https://joss.theoj.org/papers/e4e34e78e4a670f2ca9a6a97ce9d3b8e


:contributions: Please join https://github.com/cokelaer/spectrum
:contributors: https://github.com/cokelaer/spectrum/graphs/contributors
:issues: Please use https://github.com/cokelaer/spectrum/issues
:documentation: https://pyspectrum.readthedocs.io/
:Citation: Cokelaer et al, (2017), 'Spectrum': Spectral Analysis in Python, Journal of Open Source Software, 2(18), 348, doi:10.21105/joss.00348


.. image:: https://raw.githubusercontent.com/cokelaer/spectrum/main/doc/psd_all.png
    :alt: Overview of PSD methods available in Spectrum
    :width: 50%


Overview
--------

**Spectrum** contains tools to estimate Power Spectral Densities using methods
based on Fourier transform, parametric methods or eigenvalues analysis:

* **Fourier-based methods**: correlogram, periodogram and Welch estimates.
  Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are
  available (DPSS, Taylor, ...).
* **Parametric methods**: Yule-Walker, BURG, MA and ARMA, covariance and
  modified covariance methods.
* **Non-parametric (eigenanalysis) methods**: MUSIC and minimum variance
  analysis.
* **Multitapering** (MTM).

The targeted audience is diverse. Although the use of power spectrum of a
signal is fundamental in electrical engineering (e.g., radio communications,
radar), it has a wide range of applications from cosmology (e.g., detection of
gravitational waves), to music (pattern detection) or biology (mass
spectroscopy).


Quick Start
===========

The example below creates a cosine signal buried in white noise and estimates
its power spectral density using a simple periodogram:

.. code-block:: python

    from spectrum import Periodogram, data_cosine

    # generate a 1024-sample cosine at 200 Hz (amplitude 0.1) buried in white noise
    data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)

    # create the periodogram object and plot
    p = Periodogram(data, sampling=1024)
    p.plot(marker='o')

All PSD classes share the same interface: instantiate the object, run the
estimation (or let it run lazily on first access), then call ``p.plot()``.
The functional API is also available for all methods when a quick result is
needed without a full object.

See the `documentation <https://pyspectrum.readthedocs.io/>`_ for a complete
tutorial, API reference and gallery of examples.


Installation
============

**spectrum** is available on PyPI::

    pip install spectrum

and on **conda-forge**::

    conda install -c conda-forge spectrum

To install **conda** itself, see https://docs.conda.io/en/latest/miniconda.html.


Contributions
=============

Please see `GitHub <https://github.com/cokelaer/spectrum>`_ for any
issues, bugs, comments or contributions.


Changelog (summary)
====================

========== ============================================================
release    description
========== ============================================================
0.10.0     * add python 3.12 and 3.13 support
           * add sos2ss, sos2tf, etc. fix #85, #88
0.9.0      * handles new numpy API (keeping back compatibility).
           * included https://github.com/cokelaer/spectrum/pull/73
             thanks to @butala contribution to speed up fft.
           * fix rho calculation in burg algo thanks to
             https://github.com/cokelaer/spectrum/pull/82 from @cl445
           * remove warnings/deprecation related to pkgresources, numpy
             and scipy.
           * ran black through entire code.
0.8.1      * move CI to GitHub Actions
           * include Python 3.9 support
           * include PR from tikuma-lshhsc contributor to speedup
             eigenfre module
           * fix deprecated warnings
========== ============================================================


Some notebooks (external contributions)
-----------------------------------------

* http://nbviewer.ipython.org/gist/juhasch/5182528
