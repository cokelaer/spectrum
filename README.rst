SPECTRUM : Spectral Analysis in Python
==========================================

.. image:: https://badge.fury.io/py/spectrum.svg
    :target: https://pypi.python.org/pypi/spectrum

.. image:: https://secure.travis-ci.org/cokelaer/spectrum.png
    :target: http://travis-ci.org/cokelaer/spectrum

.. image:: https://coveralls.io/repos/cokelaer/spectrum/badge.png?branch=master 
    :target: https://coveralls.io/r/cokelaer/spectrum?branch=master 

.. image:: https://landscape.io/github/cokelaer/spectrum/master/landscape.png
    :target: https://landscape.io/github/cokelaer/spectrum/master

.. image:: https://anaconda.org/conda-forge/spectrum/badges/license.svg
   :target: https://anaconda.org/conda-forge/spectrum

.. image:: https://anaconda.org/conda-forge/spectrum/badges/installer/conda.svg
   :target: https://conda.anaconda.org/conda-forge

.. image:: https://anaconda.org/conda-forge/spectrum/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/spectrum


:contributions: Please join https://github.com/cokelaer/spectrum
:contributors: https://github.com/cokelaer/spectrum/graphs/contributors
:issues: Please use https://github.com/cokelaer/spectrum/issues
:documentation: http://pyspectrum.readthedocs.io/ 



.. image:: http://www.thomas-cokelaer.info/software/spectrum/html/_images/psd_all.png
    :class: align-right
    :width: 50%

**Spectrum** contains tools to estimate Power Spectral Densities using methods based on Fourier transform, Parametric methods or eigenvalues analysis:

    * The Fourier methods are based upon correlogram, periodogram and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are available (DPSS, Taylor, ...). 
    * The parametric methods are based on Yule-Walker, BURG, MA and ARMA, covariance and modified covariance methods.
    * Non-parametric methods based on eigen analysis (e.g., MUSIC) and minimum variance analysis are also implemented.
    * Multitapering is also available


The targetted audience is diverse. Although the use of power spectrum of a
signal is fundamental in electrical engineering (e.g. radio communications,
radar), it has a wide range of applications from cosmology (e.g., detection of
gravitational waves in 2016), to music (pattern detection) or biology (mass
spectroscopy).


Installation
================

::

    pip install spectrum

For Linux and MAC users, if you prefer to use conda (avoiding
 compilation of dependencies), please use::

    conda config --add channels conda-forge 
    conda install spectrum

To install the **conda** executable, please see https://www.continuum.io/downloads .

Contributions
==================

Please see `github <http://github.com/cokelaer/spectrum>`_ for any issues/bugs/comments/contributions.


Some notebooks (external contributions)
-------------------------------------------

* http://nbviewer.ipython.org/gist/juhasch/5182528
