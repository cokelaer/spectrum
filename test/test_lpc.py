from spectrum import *
from spectrum import lpc
from lpc import *
import numpy
from numpy.testing import *


def test_lpc_real():
    a, e = lpc([1,2,3,4,5,6,7,8], 3)
    assert_array_almost_equal(a, numpy.array([-0.88472690 -4.69340421e-17j,  0.01650407 +2.26160049e-17j, 0.07301958 +1.88056678e-17j]))
    assert_almost_equal(e,  9.2661256972730612)

