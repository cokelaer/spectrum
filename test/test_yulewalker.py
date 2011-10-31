from spectrum import *
from spectrum import yulewalker 
from yulewalker import *
from pylab import linspace, plot, savefig, axis, log10
from nose.tools import assert_almost_equal



#do not change. used to create figure
def test_yule():
    ar, P, c = aryule(marple_data, 15, norm='biased')
    psd = arma2psd(ar)
    return psd


# test functional checked versus octave
def test_yule_data():
    ar, v, c = aryule([1,-1,1,1,1],2, norm='biased')
    assert_almost_equal(ar[0], 0.0+0.j)
    assert_almost_equal(ar[1], -0.2+0.j)
    assert_almost_equal(v, 0.95999999999999996)
    assert_almost_equal(c[0], 0.0+0.j)
    assert_almost_equal(c[1], -0.2+0.j)

def test_pyule():
    p = pyule(marple_data, 15)
    p()
    p.plot()
    print p

#do not change. used to create figure
def create_figure():
    psd = test_yule()
    psd = cshift(psd, len(psd)/2) # switch positive and negative freq

    plot(linspace(-0.5, 0.5, 4096),
               10 * log10(psd/max(psd)))
    axis([-0.5,0.5,-60,0])
    savefig('psd_yulewalker.png')


if __name__ == "__main__":
    create_figure()
