ChangeLog Summary
===================

Version 0.7 (aug 2017)
-----------------------

* 0.7.0:

    * BUG fixes:
       * Fix https://github.com/cokelaer/spectrum/issues/38 in pburg to have the
         correct amplitude like in octave. fixed by removing the call to scale()
         function
       * Similarly all other parametric methods have been changed by adding the
         scale_by_freq argument where missing 
    * Changes:
       * remove cohere module

Version 0.6 
---------------


* 0.6.8:

    * Fix the MANIFEST

* 0.6.7:

    * refactored the requirements files (add a requirements-dev.txt) and 
      update the documentation (installation)      accordingly
    * BUG fixes: 
       * correlogram: real-data case had the data flipped
       * pmusic/pev: real-data case had the data flipped
       * fix the AKICc criteria code
    * Updates:
       * pmusic/pev: add the threshold and criteria arguments
       * more tests for the criteria and eigenfre modules
    * Changes:
       * Spectrum class: remove _correlogram method (use pcorrelogram instead)

* 0.6.6:

   * integration pull request https://github.com/cokelaer/spectrum/pull/29 from
     moritz-ritter to allow spectrum to run indepdently of matplotlib (for
     server head-less integration)

* 0.6.5:

    * minor updates to port spectrum on travis

* 0.6.4:

    * CHANGES: the bug reported in https://github.com/cokelaer/spectrum/issues/24 is
      obsolet for the reported module (pburg), which was fixed earlier but the issue
      was fixed in other module such as psd, parma, correlog
    * add LICENSE file
    * fix warning in cpp code (adding void in func() prototypes)

* 0.6.3:

    * CHANGES: portage nosetests suite to pytest
    * BUG Fixes:
    * Fix issues https://github.com/cokelaer/spectrum/issues/21 and 
      https://github.com/cokelaer/spectrum/issues/20 mostly related to
      compatibility with newest numpy version (1.12)

* 0.6.2:
    * Bug Fixes:
        * Issue #11: fixes loading mydpss library using numpy helper
        * Issue #12: Allow loading the shared library for frozen projects. Tested with py2exe.
    * Changes:
        * pmtm returns Sk_complex, weights and eigenvalues instead of just Sk

* 0.6.1:
    * BUG fixes
        * Issue #5 in pyule sampling not initialised is now fixed

* 0.6.0:
    * Code moved to github
    * plots accept the ax argument in psd module. It is a bit of a 
      hack but seems to work.

Sept 2012
----------
* 0.5.5: 
    * fix name of the libraries for mac and windows
    * change setup to manage version properly.


March 2012
--------------
* 0.5.3: add poly2lsf and lsf2poly, add tests, fix bug related to compilation of mydpss.cc
* 0.5.2: add pmtm

February 2012
--------------
* 0.5.1: add dpss wtapering windows
* 0.5.0: 
   * NPSD replaced by NFFT (qlso not correct for ARMA methods that do not have NFFT since not fourier)
   * Correlogram replaced by pcorrelogram 
   * more consistent function and class naming convention 
   * Update the entire documentation. 
* 0.4.6: fixed pylab_periodogram, documentation (installation)

January 2012
---------------

* 0.4.5: start to play with Pypi

October 2011
-----------------

* 0.4.4: Start to provide the library on the web www.assembla.com


May 2011
----------

* 0.4.3: :func:`spectrum.periodogram.pdaniell` implemented

April 2011
-----------

* 0.4.2: pcovar implemented
* 0.4.1: pmodcovar implemented
* 0.4.0: arcovar and modcovar "simplified" version. Documentation updated (tutorial, spectral_estimation, quick start...)
* 0.3.19: add linear_prediction module with codecs (eg. ac2poly, poly2rc....)
* 0.3.18 fix bug in levinson (Real data case only) and add ac2poly function.
* 0.3.17: validation of the modcovar algorithm versus the new arcovar_simplified function.
* 0.3.16: add a simplified version of arcovar called arcovar_simplified. It is 10 times faster and with a different algorithm provides the same results as arcoar, which validates the two codes!
* 0.3.15: add corrmtx function. Tested it within music algorithm
* 0.3.14: cleanup the eigen and music methods by moving the automatic order selection outside the functions.
* 0.3.13: Add AIC and MDL criteria to deal with automatic eigen values selection in pmusic and pev
* 0.3.12: test and validate the pmusic and pev pseudo spectrum.
* 0.3.11: burg and pburg  finalised
* 0.3.10: tools module cleanup and finalised
* 0.3.9:  ma fully checked and add pma validated
* 0.3.8:  minvar fully checked and add pminvar
* 0.3.7:  aryule fully checked and add pyule
* 0.3.6:  Speed up by 3 the ARMPSD (renamed to arma2psd)
* 0.3.5:  refactoring
* 0.3.4:  fix all tests and doctests
* 0.3.3:  function Daniell's periodogram implemented in module periodogram
* 0.3.2:  Create class MovingAverage, pburg, pARMA, Correlogram, Periodogram, Minvar, pma
* 0.3.1:  Cleanup MA, ARMA, BURG, MINVAR
* 0.3.0:  Create an ABC class Spectrum, a FourierSpectrum and ParametricSpectrum. 
* 0.2.4:  Finalise doc/test of the testdata module
* 0.2.3:  define a PSD class
* 0.2.2:  cleanup cholesky.py
* 0.2.1:  a new sphinx layout, 
* 0.2.0:  correlogram.py, correlation.py, levinson.py fully completed

March 2011
------------

* 31 March:
    - finalise a criteria class for AIC, FPE criteria. Incorporated it in arburg  

* 28th March:
    - First version of :func:`arcov`, :func:`aryule` and :func:`arburg` 
    - add many windows (parzen, flattop, ...).

* 22th March 2011:
    - put this doc online on thomas-cokelaer.info (fixed main links)

* 21th March 2011:
    - create psd.py defines useful class to manage Spectrum/plot
    - periodogram.py has a simple periodogram implementation equivalent to psd in pylab without overlaping. 

* 7th March 2011: 
    - add periodogram module
    - fix ARMA method in arma module
* 4th March 2011: 
    - Create first revision of spectrum package
