from spectrum import *

from pylab import linspace, log10, plot, ylim, savefig
from spectrum.periodogram import *

def test_Periodogram():
    p = Periodogram(marple_data)
    p()
    p.plot()
    print(p)

def test_periodogram():
    """check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    psd = speriodogram(marple_data)
    return psd

def test_daniell_periodogram():
    """check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    psd = DaniellPeriodogram(datasets.data_cosine(N=1024), 8, NFFT=1024)
    psd = DaniellPeriodogram(marple_data, 8, NFFT=1024)



def test_periodogram_real_vs_octave():
    # the periodogram is tested against the octave output that is "identical"
    # for the following real example
    import numpy as np
    PSDs = []
    for this in range(100):
        xx = data_two_freqs()
        p = Periodogram(xx, 4, window="hanning")
        p()
        PSDs.append(p.psd)
    M = 10*log10(np.mean(PSDs, axis=0)) 
    assert max(M)  > 10 # 10.939020375396096

    assert  np.mean(M[M < -35]) > -50
    assert  np.mean(M[M < -35]) < -40



def create_figure():
    psd = test_periodogram()
    ylim([-50,0])
    savefig('psd_periodogram.png')


if __name__ == "__main__":
    create_figure()
