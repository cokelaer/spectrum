PBURG example
================

Here is another method to estimate an AR model, based on :func:`~spectrum.burg.ARBURG` .

This example is inspired by an example found in Marple book. This is very similar
to the previous example, where you will find more explanation (see yule-Walker tutorial).



.. doctest::

    from pylab import *
    import scipy.signal
    from spectrum import *

    # Define AR filter coefficients
    a = [1, -2.2137, 2.9403, -2.1697, 0.9606];

.. doctest::

    [w,H] = scipy.signal.freqz(1, a, 256) 
    Hp = plot(w/pi, 20*log10(2*abs(H)/(2.*pi)),'r')

.. doctest::

    x = scipy.signal.lfilter([1], a, randn(256))
    AR, rho, ref = arburg(x, 4)

.. doctest::

    PSD = arma2psd(AR, rho=rho, NPSD=512)
    PSD = PSD[len(PSD):len(PSD)/2:-1]

    plot(linspace(0, 1, len(PSD)), 10*log10(abs(PSD)*2./(2.*pi)))
    xlabel('Normalized frequency (\times \pi rad/sample)')



.. plot::
    :width: 80%

    from pylab import *
    import scipy.signal
    from spectrum import *

    # Define AR filter coefficients
    a = [1, -2.2137, 2.9403, -2.1697, 0.9606];

    [w,H] = scipy.signal.freqz(1, a, 256) 
    Hp = plot(w/pi, 20*log10(2*abs(H)/(2.*pi)),'r')
    x = scipy.signal.lfilter([1], a, randn(256))
    AR, rho, ref = arburg(x, 4)
    PSD = arma2psd(AR, rho=rho, NPSD=512)
    PSD = PSD[len(PSD):len(PSD)/2:-1]
    plot(linspace(0, 1, len(PSD)), 10*log10(abs(PSD)*2./(2.*pi)))
    xlabel('Normalized frequency (\times \pi rad/sample)')
    ylabel('One-sided PSD (dB/rad/sample)')
    legend(['PSD of model output','PSD estimate of x'])


