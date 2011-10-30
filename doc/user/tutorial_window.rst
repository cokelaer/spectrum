Windowing
==========

.. contents::


In spectral analysis, it is common practice to multiply the input data by a tapering window.

Many windows are implemented and available in the :mod:`~spectrum.window` module as well as utilities
to plot the window in time and frequency domains. Some windows that have been implemented are:

.. autosummary::

    spectrum.window.window_bartlett
    spectrum.window.window_bartlett_hann
    spectrum.window.window_blackman
    spectrum.window.window_blackman_harris
    spectrum.window.window_blackman_nuttall
    spectrum.window.window_bohman
    spectrum.window.window_chebwin
    spectrum.window.window_cosine
    spectrum.window.window_flattop
    spectrum.window.window_gaussian
    spectrum.window.window_hamming
    spectrum.window.window_hann
    spectrum.window.window_kaiser
    spectrum.window.window_lanczos
    spectrum.window.window_nuttall
    spectrum.window.window_parzen
    spectrum.window.window_tukey


See :mod:`~spectrum.windows` module for a full list of windows. Note also that the :mod:`~spectrum.waveforms` provides additional waveforms/windows.

Simple window function call
----------------------------

You can explore the module to get the window function. For instance, if you look for the
Hamming window, you should find a function called :func:`~spectrum.window.window_hamming`. You 
can look at it as follows:

.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import window_hamming
    from pylab import plot

    N = 64
    w = window_hamming(N)
    plot(w)

Window Visualisation
---------------------

If you want to have a quick look at the window shape and its frequency behaviour, you can use the :func:`~spectrum.window.window_visu`:

.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import window_visu
    N = 64
    window_visu(N, 'hamming')

Window Factory
------------------


All the window are gathered within a Factory function called :func:`~spectrum.window.create_window`. The previous Hamming window can then be called using:

.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import create_window
    from pylab import plot

    N = 64
    w = create_window(N, 'hamming')
    plot(w)

The interest is that this function calls the simple function after sanity checks.


.. note:: all valid window names are stored in `spectrum.window.window_names`


Window object
--------------

Finally, there is a class :class:`~spectrum.window.Window` to further ease the manipulation of the tapering windows. 

It works as follows:


.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import Window
    N = 64
    w = Window(N, 'hamming')
    w.plot_time_freq()


You can easily access to the original data and frequency data, as well as quantities such as the equivalent noise band width:

.. doctest::

    >>> from spectrum.window import Window
    >>> N = 64
    >>> w = Window(N, 'rectangular')
    >>> w.enbw
    1.0




