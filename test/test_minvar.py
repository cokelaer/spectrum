from spectrum import *


import numpy
from spectrum import marple_data
import pylab
from numpy.testing import assert_array_almost_equal

def test_minvar_values():
    res = minvar(marple_data, 15, NFFT=16)
    assert_array_almost_equal(res[0], numpy.array([  9.34915901e-07,   9.86770648e-07,   1.38729656e-06,
         1.74456323e-06,   4.68522681e-06,   8.93086143e-05,
        -4.43764511e-06,  -3.07585044e-06,  -3.05032082e-06,
        -3.91693034e-06,  -5.80606841e-05,   6.18226003e-06,
         1.83189816e-06,   1.51951850e-06,   9.63865479e-07,
         1.00497139e-06]))

def test_minvar():
    from spectrum import tools
    res = minvar(marple_data, 15, 1.)
    psd = res[0]
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    return newpsd

def test_pminvar():
    psd = pminvar(marple_data, 15)
    psd()
    print psd
    psd = pminvar([1,2,3,4,5,6,7,-8], 2)
    psd()

def create_figure():
    res = test_minvar()
    psd = res[0]
    f = pylab.linspace(-0.5, 0.5, len(psd))
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), label='minvar 15')
    pylab.savefig('psd_minvar.png')


if __name__ == "__main__":
    create_figure()
