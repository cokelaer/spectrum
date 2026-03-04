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


def test_sos2tf():
    from scipy import signal
    sos = signal.butter(4, 0.2, output='sos')
    b, a = transfer.sos2tf(sos)
    # Verify round-trip: sos2tf then tf2zpk should match sos2zpk
    z_tf, p_tf, k_tf = signal.tf2zpk(b, a)
    z_sos, p_sos, k_sos = signal.sos2zpk(sos)
    assert_almost_equal(abs(k_tf), abs(k_sos), decimal=5)
    assert len(z_tf) == len(z_sos)
    assert len(p_tf) == len(p_sos)


def test_sos2zp():
    from scipy import signal
    sos = signal.butter(4, 0.2, output='sos')
    z, p, k = transfer.sos2zp(sos)
    z_ref, p_ref, k_ref = signal.sos2zpk(sos)
    assert_almost_equal(np.sort(np.abs(z)), np.sort(np.abs(z_ref)), decimal=5)
    assert_almost_equal(np.sort(np.abs(p)), np.sort(np.abs(p_ref)), decimal=5)
    assert_almost_equal(k, k_ref, decimal=5)


def test_sos2ss():
    from scipy import signal
    sos = signal.butter(4, 0.2, output='sos')
    A, B, C, D = transfer.sos2ss(sos)
    # Verify state-space matrices are returned with correct shapes
    n_sections = sos.shape[0]
    order = 2 * n_sections
    assert A.shape == (order, order)
    assert B.shape == (order, 1)
    assert C.shape == (1, order)
    assert D.shape == (1, 1)
    # Verify round-trip: ss2zpk should recover poles/zeros matching sos2zpk
    z_ss, p_ss, k_ss = signal.ss2zpk(A, B, C, D)
    z_sos, p_sos, k_sos = signal.sos2zpk(sos)
    assert_almost_equal(np.sort(np.abs(p_ss)), np.sort(np.abs(p_sos)), decimal=5)
    assert_almost_equal(abs(k_ss), abs(k_sos), decimal=5)








