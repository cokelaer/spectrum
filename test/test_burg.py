from spectrum.arma import arma2psd
from spectrum.burg import *
from spectrum.datasets import marple_data, data_cosine
from spectrum.tools import cshift
from spectrum.criteria import *

from numpy.testing import *
import numpy
import pylab


def test_arburg_functional():
    ar, P, k = arburg(marple_data, order=15)
    PSD = arma2psd(ar)
    newpsd = cshift(PSD, len(PSD)/2) # switch positive and negative freq
    return newpsd


def test_arburg_real_output():
    a,b,c = arburg([1,2,3,4,5,6,7,-8], 4)
    assert_array_almost_equal(a, [-3.2581e-01,   3.4571e-04,   3.3790e-02,   9.8853e-02], decimal=5)
    assert_almost_equal(b, 22.44687609887804)
    assert_array_almost_equal(c, [-0.326531,   0.022568, 0.066648,  0.098853], decimal=5)


def test_arburg_imag_output():
    a, b, c = arburg(marple_data, 15)


    a_e, b_e, c_e = (numpy.array([ 2.70936368 -0.77610302j,  5.17482864 -2.73293024j,
        7.03527787 -6.15070038j,  7.89423853-10.20591369j,
        6.84853701-14.07469247j,  4.56915619-16.84486008j,
        1.32687590-18.13284671j, -1.87811360-17.49937286j,
       -4.64976221-15.05888331j, -6.22557823-11.25070227j,
       -6.28367510 -6.93498375j, -4.89652279 -3.24910899j,
       -2.99758653 -0.8736847j , -1.32183647 +0.04527281j,
       -0.35565856 +0.14754881j]),
     0.0054379699760549929,
     numpy.array([-0.18570222-0.87179346j,  0.26402371-0.5190592j ,
        0.07162311-0.46372011j,  0.44463099+0.05080174j,
       -0.02634972-0.14691215j,  0.19255061-0.37032848j,
       -0.25994598-0.55924338j, -0.20237974-0.23641516j,
       -0.40546748-0.40598876j, -0.47824854-0.42553068j,
       -0.51507096-0.49435948j, -0.32530245-0.49134098j,
       -0.21950049-0.37261937j, -0.28613904-0.0921211j ,
       -0.35565856+0.14754881j]))


    assert_array_almost_equal(a, a_e)
    assert_almost_equal(b,  b_e)
    assert_array_almost_equal(c, c_e)



def test_burg_criteria():
    ar, P, k = arburg(marple_data, order=15, criteria='AIC')
    ar, P, k = arburg(marple_data, order=15, criteria='AICc')
    ar, P, k = arburg(marple_data, order=15, criteria='KIC')
    ar, P, k = arburg(marple_data, order=15, criteria='MDL')
    ar, P, k = arburg(marple_data, order=15, criteria='FPE')
    ar, P, k = arburg(marple_data, order=15, criteria='AKICc')

def test_pburg():
    p = pburg(marple_data, 15, NFFT=4096)
    p()
    p.plot(sides='centerdc')
    print p
    # test real case
    p = pburg(data_cosine(), 15, NFFT=4096)
    p()

def create_figure():
    psd = test_burg()
    pylab.plot(pylab.linspace(-0.5, 0.5, len(psd)),
               10 * pylab.log10(psd/max(psd)))
    pylab.axis([-0.5,0.5,-60,0])
    pylab.savefig('psd_burg.png')


if __name__ == "__main__":
    create_figure()
