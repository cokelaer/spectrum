Windowing
==========

.. contents::


In spectral analysis, it is common practice to multiply the input data by a tapering window.

Many windows are implemented and available in the :mod:`~spectrum.window` module as well as utilities
to plot the window in time and frequency domains. Some windows that have been implemented are:

.. autosummary::

    spectrum.window.window_bartlett
    spectrum.window.window_blackman
    spectrum.window.window_gaussian
    spectrum.window.window_hamming
    spectrum.window.window_hann
    spectrum.window.window_kaiser
    spectrum.window.window_lanczos
    spectrum.window.window_nuttall
    spectrum.window.window_tukey



See :mod:`~spectrum.window` module for a full list of windows. Note also that the :mod:`~spectrum.waveform` provides additional waveforms/windows.


Window object
--------------

There is a class :class:`~spectrum.window.Window` that ease the manipulation of the tapering windows. It works as follows:


.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import Window
    N = 64
    w = Window(N, 'hamming')
    w.plot_time_freq()

where `N` is the length of the desired window, and "hamming" is the name. There are a lot of different windows, some of them require arguments. For example, the blackman window require an `alpha` argument::

    w = Window(64, 'blackman', alpha=1)

From the object, you can easily access to the window data (`w.data`) and frequency (`w.frequencies`), as well as quantities such as the equivalent noise band width:

.. doctest::

    >>> from spectrum.window import Window
    >>> N = 64
    >>> w = Window(N, 'rectangular')
    >>> w.enbw
    1.0

To have a list of valid names, omit the name. It should raise an error with the list of valid names. Alternatively, type::

    window_names.keys()

Finally, when a window require arguments, you need to know their names (e.g., in the blackman example above, the `alpha` parameter is required).

The only way to get this information is to look at the function `window_<name>` (e.g. window_blackman) and type::

    window_blackman?



Simple window function call
----------------------------

You can explore the module to get the window function. For instance, if you look for the Hamming window, you should find a function called :func:`~spectrum.window.window_hamming`. You can look at it as follows:

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

If you do not want the object approach, you may want to use the Factory function called :func:`~spectrum.window.create_window` (the :class:`Window` class relies on this function). The previous Hamming window can be called using:

.. plot::
    :width: 80%
    :include-source:

    from spectrum.window import create_window
    from pylab import plot

    N = 64
    w = create_window(N, 'hamming')
    plot(w)

