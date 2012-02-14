All PSD methods
================

This example is used to generate the front image. It shows how to use the different PSD classes that can be found in **Spectrum**.

.. plot::
    :include-source:
    :width: 80%


    from spectrum import *
    from pylab import legend, ylim
    data = marple_data
    norm = True
    sides = 'centerdc'

    # MA method
    p = pma(marple_data, 15, 30, NFFT=4096)
    p(); p.plot(label='MA (15, 30)', norm=norm, sides=sides)

    # ARMA method
    p = parma(marple_data, 15, 15, 30, NFFT=4096)
    p(); p.plot(label='ARMA(15,15)', norm=norm, sides=sides)

    # yulewalker
    p = pyule(data, 15, norm='biased', NFFT=4096)
    p(); p.plot(label='YuleWalker(15)', norm=norm, sides=sides)

    #burg method
    p = pburg(data, order=15, NFFT=4096)
    p(); p.plot(label='Burg(15)', norm=norm, sides=sides)

    #covar method
    p = pcovar(data, 15, NFFT=4096)
    p(); p.plot(label='Covar(15)', norm=norm, sides=sides)

    #modcovar method
    p = pmodcovar(data, 15, NFFT=4096)
    p(); p.plot(label='Modcovar(15)', norm=norm, sides=sides)

    # correlagram
    p = pcorrelogram(data, lag=15, NFFT=4096)
    p(); p.plot(label='Correlogram(15)', norm=norm, sides=sides)

    #minvar
    p = pminvar(data, 15, NFFT=4096)
    p(); p.plot(label='minvar (15)', norm=norm, sides=sides)

    #music
    p = pmusic(data, 15, 11, NFFT=4096)
    p(); p.plot(label='music (15, 11)', norm=norm, sides=sides)

    #ev
    p = pev(data, 15, 11, NFFT=4096)
    p(); p.plot(label='ev (15, 11)', norm=norm, sides=sides)


    legend(loc='upper left', prop={'size':10}, ncol=2)
    ylim([-80,10])

