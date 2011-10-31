from spectrum import transfer
from transfer import *
from nose.tools import assert_almost_equal

def test_tf2zpk():
    from scipy import signal
    [b,a] = signal.butter(3,.4);
    z,p,k = transfer.tf2zpk(b,a)
    assert_almost_equal(k, 0.09853, 4)
    assert_almost_equal(z[0], -1.00000293 +5.07009147e-06j, 4)
    assert_almost_equal(z[1], -1.00000293 -5.07009147e-06j, 4)
    assert_almost_equal(z[2], -0.99999415 +0.00000000e+00j, 4)
    assert_almost_equal(p[0],  0.20942804+0.55819948j)
    assert_almost_equal(p[1], 0.20942804-0.55819948j)
    assert_almost_equal(p[2], 0.15838444+0.j)


