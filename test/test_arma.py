from spectrum import *
from spectrum import arma
from arma import *
from pylab import linspace, log10, plot, ylim, savefig
from nose.tools import assert_almost_equal
from spectrum.arma import arma_estimate, arma2psd, parma, pma, ma
from numpy.testing import assert_array_almost_equal, assert_almost_equal
import numpy

def test_arma_values():
    a, b, rho = arma_estimate(marple_data, 15, 15, 30)
    assert_almost_equal(rho, 0.20050144053393698, decimal=6)

    assert_array_almost_equal(a, numpy.array([ 
        1.47857824-0.16358208j,  4.32139091-0.86231938j,
        6.04115773-1.77030183j,  6.09285854-3.96367752j,
        4.70699008-3.27199141j,  3.45467782-1.59183506j,
        3.11230094-1.06510595j,  1.55237009+1.09800024j,
        1.05148353+2.2720917j ,  1.68042547+4.9737292j ,
        3.22899406+6.39981425j,  3.16557650+5.92783737j,
        3.47120865+5.48246963j,  2.79508215+3.3238971j ,
        2.13174602+1.51034329j  
        ]), decimal=4)
    

def test_arma():
    """arma, check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    a,b, rho = arma_estimate(marple_data, 20,20,40)
    psd = arma2psd(A=a,B=b, rho=rho)
    return psd

def create_figure_arma():
    psd = test_arma()
    psd = cshift(psd, len(psd)/2) # switch positive and negative freq
    plot(linspace(-0.5, 0.5, 4096), 10 * log10(psd/max(psd)))
    ylim([-50,0])
    savefig('psd_arma.png')

def create_figure_ma():
    psd = test_ma()
    psd = cshift(psd, len(psd)/2) # switch positive and negative freq
    plot(linspace(-0.5, 0.5, 4096), 10 * log10(psd/max(psd)))
    ylim([-50,0])
    savefig('psd_ma.png')

def test_ma():
    """ma PSD. check that rho is correct (appendix 10.A )and reproduce figure 10.2"""
    b, rho = ma(marple_data, 15, 30)
    assert_almost_equal(rho, 0.21432, decimal=4)

    assert_almost_equal(b[0], -0.25150803+0.67246418j, decimal=6)
    assert_almost_equal(b[1], -0.68612023+0.14571702j, decimal=6)

    """-0.25150803+0.67246418j, -0.68612023+0.14571702j,
        0.02061484-0.52246411j,  0.11444091-0.19157961j,
        0.36592370+0.09885689j,  0.00556917+0.4330789j ,
       -0.40634639+0.04854752j, -0.10092740-0.42813962j,
        0.26013726-0.16101382j,  0.25119793+0.16711825j,
       -0.13448885+0.21202256j, -0.16125290+0.10804393j,
       -0.03402254-0.18015694j,  0.08780647-0.1146388j ,
        0.02294750+0.08411391j
    """
    psd = arma2psd(B=b, rho=rho)
    return psd


def test_arma2psd():
    psd = arma2psd([0.5], NPSD=16, norm=True)*4
    assert_array_almost_equal(psd, numpy.array([ 0.44444444,  0.46000709,  0.51095832,  0.61248861,  0.8       ,
            1.15298155,  1.84198285,  3.06635155,  4.        ,  3.06635155,
            1.84198285,  1.15298155,  0.8       ,  0.61248861,  0.51095832,
            0.46000709]))



def test_parma():
    p = parma(marple_data, 4, 4, 30, NFFT=4096)
    p()
    p.plot()
    print p

def test_moving_average_class():
    p = pma(marple_data, 15, 30, NFFT=4096)
    p()
    print p
    p = pma(data_cosine(N=1024), 15, 30, NFFT=4096)

def create_figure_ma():
    psd = test_ma()
    psd = cshift(psd, len(psd)/2) # switch positive and negative freq
    plot(linspace(-0.5, 0.5, 4096), 10 * log10(psd/max(psd)))
    ylim([-50,0])
    savefig('psd_ma.png')


if __name__ == "__main__":
    create_figure_ma()
    create_figure_arma()
