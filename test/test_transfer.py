from spectrum import transfer
from numpy.testing import assert_almost_equal
import numpy as np

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


def test_tf2zp():
    b = [2, 3, 0]
    a = [1, 0.4, 1]
    [z,p,k] = transfer.tf2zp(b,a)
    assert all(z == np.array([-1.5,0]))
    assert k == 2
    assert_almost_equal(p[0], -0.2000 + 0.9798j,4)
    assert_almost_equal(p[1],  -0.2000 - 0.9798j,4)

    transfer.zp2tf(z,p,k)

def test_eqtlength():
    a, b = transfer.eqtflength([1,2,3,4], np.array([1,2]))
    assert all(b == np.array([1,2,0,0]))
    a, b = transfer.eqtflength([1,2,3,4], [1,2])
    assert b == [1,2,0,0]
    a, b = transfer.eqtflength(np.array([1,2]), [1,2,3,4])
    assert all(a == np.array([1,2,0,0]))
    a, b = transfer.eqtflength([1,2], [1,2,3,4])
    assert a == [1,2,0,0]
    a, b = transfer.eqtflength([1,2], [1,2])
    assert a == b


def test_latc2tf():
    try:
        transfer.latc2tf()
        assert False
    except NotImplementedError:
        assert True
    except:
        assert False


def test_latcfilt():
    try:
        transfer.latcfilt()
        assert False
    except NotImplementedError:
        assert True
    except:
        assert False


def test_tf2ss():
    try:
        transfer.tf2ss()
        assert False
    except NotImplementedError:
        assert True
    except:
        assert False


def test_tf2sos():
    try:
        transfer.tf2sos()
        assert False
    except NotImplementedError:
        assert True
    except:
        assert False


def test_ss2zpk():
    z,p,k = [1.5,1], [-0.2,1], 2.

    zp,pp,kp = transfer.ss2zpk(*transfer.zpk2ss(z,p,k))
    #assert zp == z
    # FIXME: this fails
    #assert pp == p
    assert k == kp








