What is the Spectrum object ?
===============================

.. module:: spectrum.psd


Normally Users should not be bother by the classes used. For instance if you use the pburg class to compute a PSD estimate base on the Burg method, you just nee to use :class:`~spectrum.burg.pburg`. Indeed, the normal usage to estimate a PSD is to use the PSD estimate starting with the letter `p` such as parma, pminvar, pburg, (exception: use Periodogram instead of pPeriodogram). 


Yet, it may be useful for some advanced users and developers to know that all PSD estimates are based upon the :class:`Spectrum` class (used by specialised classes such as :class:`FourierSpectrum` and :class:`ParametricSpectrum`).

The following example shows how to use :class:`Spectrum`. First, let us create a Spectrum instance (first argument is the time series/data)::

    from spectrum import *
    p = Spectrum(data_cosine(), sampling=1024)

Some information are stored and can be retrieved later on::

    p.N
    p.sampling


However, for now it contains no information about the PSD estimation method. For instance, if you type::

    p.psd

it should return a warning message telling you that the PSD has not yet been computed. You can compute it either independantly, and set the `psd` attribute manually::

    psd = speriodogram(p.data)

or you can associate a function to the `method` attribute::

    p.method = minvar

and then call the function with the proper optional arguments::

    p(15, NFFT=4096)

In both cases, the PSD is now saved in the `psd` attribute.

Of course, if you already know the method you want to use, then it is much simpler to call the appropriate class directly as shown in previous sections and examples::

    p =  pminvar(data_cosine(), 15)
    p()
    p.plot()


.. plot:: 
    :width: 80%

    from spectrum import *
    p =  pminvar(data_cosine(), 15)
    p()
    p.plot()

