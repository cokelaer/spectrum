#from spectrum.linear_prediction import *
import numpy
from numpy.testing import *

ri = [5.0000, -1.5450, -3.9547, 3.9331, 1.4681, -4.7500]
ki = [0.3090 ,   0.9800    ,0.0031  ,  0.0082 ,  -0.0082];
ai = [1, 0.6147,  0.9898 , 0.0004,  0.0034, -0.0077]

def _test_ac2poly():
    [a,efinal] = ac2poly(ri)
    assert_array_almost_equal(a, ai, decimal=4)
    assert_almost_equal(efinal, 0.1791, decimal=4)

def _test_rc2poly():
    a, e = rc2poly(ki, 1)
    assert_array_almost_equal(a , numpy.array([1.00000000e+00,   6.14816180e-01,   9.89881431e-01,
         2.42604054e-05,   3.15795596e-03,  -8.20000000e-03]), decimal=4)
    assert_almost_equal(e, 0.035813791455383194)

def _test_rc2ac():
    R = rc2ac(ki, 1)
    assert_array_almost_equal(R, numpy.array([ 1.00000000+0.j, -0.30900000+0.j, -0.79094762-0.j,  0.78662653-0.j, 0.29362938-0.j, -0.95000010-0.j]))


def _test_ac2rc():
    a, r0 = ac2rc(ri)
    assert_almost_equal(r0, 5)


def _test_poly2ac():
    a, efinal = ac2poly(ri)
    r = poly2ac(a, efinal)
    assert_array_almost_equal(a, ai, decimal=4)


def _test_poly2rc():
    a, efinal = ac2poly(ri)
    k = poly2rc(numpy.insert(a, 0, 1), efinal)
    #assert_array_almost_equal(ki, k)
    #\numpy.array([ 0.309     ,  0.97999158,  0.00302085,  0.00818465, -0.00770967]))



def _test_lar2rc():
    lar2rc(.1)
    lar2rc([.1,.2])

def _test_rc2lar():
    #rc2lar(0.5)
    rc2lar([0.5,0.9])
    try:
        rc2lar([2])
        assert False
    except:
        assert True

def _test_is2rc():
    is2rc([1,2])

def _test_rc2is():
    #rc2is(0.5)
    rc2is([0.5,0.9])
    try:
        rc2is([2])
        assert False
    except:
        assert True
