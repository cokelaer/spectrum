from spectrum import *
from spectrum import modcovar 
from modcovar import *
import pylab
from numpy.testing import *

def test_modcovar_ip_null():
    a, p, pv = modcovar_marple(marple_data, 0)
    assert_almost_equal(p, 1.780459894489305)
    assert len(a) == 0
    assert len(pv) == 0


def test_modcovar():
    a, p, pv = modcovar_marple(marple_data, 15)
    PSD = arma2psd(a, sides='centerdc')
    return PSD


def test_covar_simplified():
    a, b, c = modcovar_marple(marple_data, 15)
    a2, e2 = modcovar(marple_data, 15)
    assert_array_almost_equal(a[0:15], a2)  # af contains zeros after order=15

def test_pmodcovar():
    # test real data
    p  = pmodcovar(data_cosine(), 15, scale_by_freq=True)
    p()
    print p

    # and complex data
    p  = pmodcovar(marple_data, 15, NFFT=4096)
    p()

    # return psd for the create_figure function
    p.sides = 'centerdc'
    return p.psd

def create_figure():
    newpsd = test_modcovar()

    pylab.plot(pylab.linspace(-0.5, 0.5, 4096),
               10 * pylab.log10(newpsd/max(newpsd)))
    pylab.axis([-0.5,0.5,-60,0])
    pylab.savefig('psd_modcovar.png')
 
if __name__ == "__main__":
    create_figure()
