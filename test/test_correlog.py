from spectrum import *

from nose.tools import assert_almost_equal
from pylab import *
from numpy.testing import assert_array_almost_equal

def test_correlog():
    psd = CORRELOGRAMPSD(marple_data, marple_data, lag=15)
    assert_almost_equal(psd[0],      0.138216970)
    assert_almost_equal(psd[1000-1], 7.900110787)
    assert_almost_equal(psd[2000-1], 0.110103858)
    assert_almost_equal(psd[3000-1], 0.222184134)
    assert_almost_equal(psd[4000-1], -0.036255277)
    assert_almost_equal(psd[4096-1], 0.1391839711)
    return psd

def test_correlog_auto_cross():
    """Same as test_correlog but x and y provided"""
    psd1 = CORRELOGRAMPSD(marple_data, lag=16)
    psd2 = CORRELOGRAMPSD(marple_data, marple_data, lag=16)
    assert_array_almost_equal(psd1, psd2)

    psd1 = CORRELOGRAMPSD(marple_data, lag=16, correlation_method='CORRELATION')
    psd2 = CORRELOGRAMPSD(marple_data, marple_data, lag=16, correlation_method='CORRELATION')
    assert_array_almost_equal(psd1, psd2)

def test_correlog_correlation_method():
    """test correlogramPSD playing with method argument"""
    psd1 = CORRELOGRAMPSD(marple_data, lag=16, correlation_method='CORRELATION')
    psd2 = CORRELOGRAMPSD(marple_data, marple_data, lag=16, correlation_method='xcorr')
    assert_array_almost_equal(psd1, psd2)

def test_pcorrelogram_class():
    p = pcorrelogram(marple_data, lag=16)
    p()
    print p

def create_figure():
    psd = test_correlog()
    f = linspace(-0.5, 0.5, len(psd))

    psd = cshift(psd, len(psd)/2)
    plot(f, 10*log10(psd/max(psd)))
    savefig('psd_corr.png')

if __name__ == "__main__":
    create_figure()


