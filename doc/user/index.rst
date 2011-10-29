User Guide
===========

QuickStart (Periodogram example)
--------------------------------

The simplest way to start with **Spectrum** is to import everything from the library:

.. doctest::

    from spectrum import *

**Spectrum** provides some function that creates data samples for you. For instance, the :func:`~spectrum.datasets.data_cosine` function creates a noisy cosine. It works as follows::

    data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)

`data` contains a cosine signal with a frequency of 200Hz and buried in white noise (amplitude 0.1). The data has a length of N and the sampling is 1024Hz.

Now, we can analyse this data using one of the Power Spectrum Estimation method. You can either use a functional or object approach. For now, we will use the object approach because it provides more robustness and additional tools (e.g., plotting). So, let us create a simple periodogram::

    p = Periodogram(data, sampling=1024) 
    p() # now you run the estimation
    p.plot(marker='o')  # standard matplotlib options are accepted


.. plot::
    :width: 80%

    from spectrum import *
    data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)
    p = Periodogram(data, sampling=1024) #here you just define the PSD estimate 
    p() # now you run the estimation
    p.plot(marker='o')

Since the data is purely real, the PSD (stored in p.psd) is a onesided PSD, with positive frequencies only. If the data were complex, the two-sided PSD would have been computed and plotted. For the real case, you can still plot a two-sided PSD by setting the sides option manually::

    p.plot(sides='twosided')

You can also look at a centered PSD around the zero frequency::

    p.plot(sides='centerdc')

.. warning:: By convention, the :attr:`psd` attribute contain the default PSD (either one-sided for real data or two-sided for the complex data).

Since **p** is an instance of Periodogram, you can introspect the object to obtain diverse information such as the original data, the sampling, the PSD itself and so on::

   p.psd # contains the PSD values
   p.frequencies() returns a list of frequencies
   print(p) # prints some information about the PSD.


The object approach versus functional approach (ARMA example)
--------------------------------------------------------------

Object approach
~~~~~~~~~~~~~~~~~~
In the previous section, we've already illustrated the object approach using a Fourier-based method with the simple periodogram method. In addition to the Fourier-based PSD estimates, **Spectrum** also provides parametric-based estimates. Let us use :func:`~spectrum.arma.parma` class as a second illustrative example of the object approach:

.. doctest::

    from spectrum import parma


Many functionalities available in **Spectrum** are inspired by methods found in [Marple]_. The data sample used in most of the examples in this reference are available and can be imported as follows (a 64 complex data samples)::

    from spectrum import marple_data

The class :class:`~spectrum.arma.parma` allows to create an ARMA model and to plot the PSD, similarly to the previous example (Periodogram). First, we need to create the object::

    p = parma(marple_data, 15, 15, 30, NFFT=4096)

where 15,15 and 30 are arguments of the ARAM model (see :class:`spectrum.parma`).

Then, computation and plot can be performed::

    p()
    p.plot(norm=True, color='red', linewidth=2)

.. plot::
    :width: 80%

    from spectrum import parma, marple_data
    p = parma(marple_data, 15, 15, 30, NFFT=4096)
    p() # now you run the estimation
    p.plot(norm=True, color='red', linewidth=2) # same options as pylab.plot

Since the data is complex, the PSD (stored in p.psd) is a twosided PSD. Note also that all optional arguments accepted by matplotlib function are also available in this implementation. 


Functional approach
~~~~~~~~~~~~~~~~~~~~
The object-oriented approach can be replaced by a functional one if required. Nevertheless, as mentionned earlier, this approach required more expertise and could easily lead to errors. Here is an example that is identical to the previous piece of code. First, we need two functions, one for the estimation, one for the PSD computation (and plotting)::

    from spectrum.arma import arma_estimate, arma2psd

In order to extract the autoregressive coefficients (AR) and Moving average coefficients (MA), the :func:`~spectrum.arma.arma_estimate` can be used::

    ar, ma, rho = arma_estimate(marple_data, 15, 15, 30)

When you have AR and/or MA parameters, the :func:`~spectrum.arma.arma2psd` function creates a two-sided PSD for you. Combining the results with plotting routines from Pylab, we the following scripts creates the plot:

.. plot::
    :include-source:
    :width: 80%

    from spectrum import arma_estimate, arma2psd, marple_data
    from pylab import *
    ar, ma, rho = arma_estimate(marple_data, 15, 15, 30)
    psd = arma2psd(ar, ma, rho=rho, NPSD=4096)
    plot(10*log10(psd/max(psd)))
    axis([0, 4096, -80, 0])
    xlabel('Frequency')
    ylabel('power (dB)')
    grid(True)

.. note::

    #. The parameter 30 is the correlation lag that should be twice as much as the required AR and MA coefficient number (see reference guide for details). 
    #. Here we plot the PSD manually, and normalise it so as to use dB units (10*log10)
    #. Since the data are complex data, the default plot is a two-sided PSD.
    #. The frequency vector is not provided. 

