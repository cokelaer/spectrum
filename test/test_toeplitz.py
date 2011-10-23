import numpy
from nose.tools import assert_almost_equal
from spectrum import toeplitz
from spectrum.toeplitz import * #so that nosetests scans te docstrings as well
from spectrum.toeplitz import TOEPLITZ

def test_toeplitz():
    T0 = 3+0j
    TC = numpy.array([-2.+.5j, .7-1j ])
    TR = numpy.array([-.2-.4j, .3-.6j])
    Z = numpy.array([1.+3.j, 2.-1.j, .5+.8j])
    X = TOEPLITZ(T0, TC, TR, Z)
    return X

    expected_X = numpy.array([ 0.23518549+1.24372203j,  1.03025807+0.53575117j,  0.47334662+0.24031779j])
    assert_almost_equal(X.all(), expected_X.all())


def test_hermtoep():
    T0=3.
    T = numpy.array([-2.+.5j, .7-1j ])
    Z = numpy.array([1.+3.j, 2.-1.j, .5+.8j])
    X = HERMTOEP(T0, T, Z)
    return X
    expected_X = numpy.array([ 2.29697452+1.69904459j,  3.31847134+1.29697452j, 1.49283439+0.94745223j])
    assert_almost_equal(X.all(), expected_X.all())

