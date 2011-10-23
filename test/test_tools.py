from numpy.testing import assert_almost_equal, assert_array_almost_equal
from numpy import *
from spectrum.tools import *


def test_nextpow2():
    res = nextpow2([255,256,257])
    assert res[0] == 8.
    assert res[1] == 8.
    assert res[2] == 9.

def test_pow2db_db2pow():
    p1 = 10.
    x = pow2db(p1)
    p2 = db2pow(x)
    assert_almost_equal(p1, p2)

def test_mag2db_db2mag():
    p1 = 10.
    x = mag2db(p1)
    p2 = db2mag(x)
    assert_almost_equal(p1, p2)

def test_cshift():
    a = [1,2,3,4]
    b = cshift(a, 2)
    assert_array_almost_equal([3,4,1,2], b)
   

def _test_twosided_zerolag():
    data = [3,2,1]
    zerolag = 4
    res = twosided_zerolag(data, zerolag)
    assert_array_almost_equal(array([1, 2, 3, 4, 3, 2, 1]), res)

def test_twosided():
    a = [1,2,3]
    b = twosided(a)
    assert_array_almost_equal(b, array([3, 2, 1, 2, 3]))

def _test_swap_sides():
    x = [-2, -1, 1, 2]
    b = swapsides(x)
    assert_array_almost_equal(b, array([2, -2, -1]))

def _test_fftshift():
    x = [1,2,3,4,5]
    y = fftshift([1,2,3,5,4])
    assert_array_almost_equal(y, array([5, 4, 1, 2, 3]))


def test_onesided_twosided():
    x = [10, 2, 3, 4, 6]
    y = onesided_2_twosided(x)
    x2 = twosided_2_onesided(y)
    assert_array_almost_equal(x, x2)
    
def test_centeddc_twosided():
    x = [10, 2, 3, 4, 4, 3, 2, 20]
    y = centerdc_2_twosided(x)
    x2 = twosided_2_centerdc(y)
    assert_array_almost_equal(x, x2)
