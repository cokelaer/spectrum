from spectrum import *

from pylab import linspace, log10, plot, ylim, savefig
from nose.tools import assert_almost_equal
from spectrum.periodogram import *


def test_Periodogram():
    p = Periodogram(marple_data)
    p()
    p.plot()
    print p

def test_periodogram():
    """check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    psd = speriodogram(marple_data)
    return psd

def test_daniell_periodogram():
    """check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    psd = DaniellPeriodogram(datasets.data_cosine(N=1024), 8, NFFT=1024)
    psd = DaniellPeriodogram(marple_data, 8, NFFT=1024)


def create_figure():
    psd = test_periodogram()
    ylim([-50,0])
    savefig('psd_periodogram.png')


if __name__ == "__main__":
    create_figure()
