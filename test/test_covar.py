from spectrum import *
import pylab
from nose.tools import assert_almost_equal
from numpy.testing import *

def test_covar_null_ip():
    af, pf, ab, pb, c = arcovar_marple(marple_data, 0)
    assert_almost_equal(pf, 1.7804598944893046)

def test_covar_15_ip():
    af, pf, ab, pb, pbv = arcovar_marple(marple_data, 15)
    assert_almost_equal(pf, 0.0031358526195905032)
    assert_almost_equal(pb, 0.0026095580050847235)
    assert_array_almost_equal(af[0:15], array([  3.14064291e+00 -0.53085796j,   6.71499124e+00 -2.02047795j,
         1.06218919e+01 -4.91215366j,   1.40604378e+01 -8.88144555j,
         1.56600743e+01-13.2925649j ,   1.52808636e+01-17.26357445j,
         1.29553371e+01-20.19441487j,   9.56479043e+00-21.35967801j,
         5.76086019e+00-20.39407074j,   2.35478080e+00-17.25236853j,
        -1.39883911e-02-12.63099132j,  -1.01307484e+00 -7.71542788j,
        -1.00735874e+00 -3.71449987j,  -5.47782956e-01 -1.24481265j,
        -1.63739470e-01 -0.22820697j]))
    assert_array_almost_equal(ab[0:15], array([  3.06854326 +0.4396126j ,   6.52836187 +1.85223579j,
        10.14250939 +4.53484335j,  13.27104933 +8.16295648j,
        14.65282324+12.10370542j,  14.30283278+15.67072521j,
        12.13984749+18.32533332j,   9.02885933+19.34952244j,
         5.49933445+18.38815454j,   2.39313549+15.41172794j,
         0.23240843+11.16952573j,  -0.69430878 +6.74812076j,
        -0.75349882 +3.21552564j,  -0.42710881 +1.07407686j,
        -0.13625884 +0.18990667j]))
    assert_array_almost_equal(pbv, array([23.002882564886164,
         14.963158025030376,     11.46060060362683,
         8.8047876198403294,     8.464718707735825,
         6.7595928955003961,     3.9194229830412644,
         3.4283223276191257,     2.2528330561384045,
         1.174361182536527,     0.53260425403862111,
         0.30138304540853789,     0.1893577453852136,
         0.13685257356088598]))


def test_covar_simplified():
    af, pf, ab, pb, pv = arcovar_marple(marple_data, 15)
    a2, e2 = arcovar(marple_data, 15)
    assert_array_almost_equal(af[0:15], a2)  # af contains zeros after order=15


def test_covar():
    af, pf, ab, pb, pv = arcovar_marple(marple_data, 15)
    PSD = arma2psd(af)

    newpsd = cshift(PSD, len(PSD)/2) # switch positive and negative freq
    return newpsd

def test_pcovar():
    p = pcovar(data_cosine(), 15, NFFT=4096, scale_by_freq=True)
    p()
    p = pcovar(marple_data, 15, NFFT=4096)
    p()
    print p.get_converted_psd('centerdc')
    return p.psd


def create_figure():
    psd = test_pcovar()
    pylab.axis([-0.5,0.5,-60,0])
    pylab.savefig('psd_covar.png')


if __name__ == "__main__":
    create_figure()
