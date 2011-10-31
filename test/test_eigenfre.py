from spectrum import *
from spectrum import eigenfre
from eigenfre import *

from pylab import plot, linspace, log10, savefig
import numpy
from numpy.testing import *



def test_pmusic():
    p = pmusic(marple_data, 15, NSIG=11)    
    p()
    p = pmusic(data_cosine(), 15, NSIG=11)    
    p()
    print p


def test_eigenfre_music():

    psd, s = music(marple_data, 15, NSIG=11)
    return psd


def test_eigenfre_ev():
    psd, s = ev(marple_data, 15, NSIG=11)
    assert_array_almost_equal(s, numpy.array([  4.45510959e+01,   1.01451096e+01,   8.37309134e+00,
         7.17637043e+00,   6.62545637e+00,   5.83043837e+00,
         4.16284271e+00,   3.69224764e+00,   3.64345761e+00,
         3.07519938e+00,   2.05618798e+00,   1.53143913e+00,
         8.21242005e-01,   1.10463229e-01,   1.02225490e-02]))
    return psd

def test_eigen_parameters():
    psd, s = ev(data_cosine(), 15)
    psd, s = ev(data_cosine(), 15, NSIG=11)
    psd, s = ev(data_cosine(), 15, threshold=2)

def create_figure():
    psd = test_eigenfre_music()
    f = linspace(-0.5, 0.5, len(psd))
    plot(f, 10 * log10(psd/max(psd)), '--',label='MUSIC 15')
    savefig('psd_eigenfre_music.png')

    psd = test_eigenfre_ev()
    f = linspace(-0.5, 0.5, len(psd))
    plot(f, 10 * log10(psd/max(psd)), '--',label='EV 15')
    savefig('psd_eigenfre_ev.png')


if __name__ == "__main__":
    create_figure() 
