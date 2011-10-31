
from spectrum import *
from spectrum import cholesky
from cholesky import *
from numpy.testing import assert_array_almost_equal
import numpy

#TOEPLITZ
#example marple app 3.A
A  = numpy.matrix([[2+0.j, .5-0.5j,-.2+.1j],[.5+.5j,1,.3-0.2j],[-.2-.1j,.3+.2j,.5]], dtype=complex)

#a_cholesky = numpy.array([2+0.j, .5-0.5j, 1., -.2+.1j,.3-0.2j,.5], dtype=complex)
B = numpy.array([1+3j,2-1j,.5+.8j], dtype=complex)
#should return 
sol_cholesky = numpy.array([ 0.95945946+5.25675676j, 4.41891892-7.04054054j, 
-5.13513514+6.35135135j])


def test_cholesky():
    a = CHOLESKY(A, B, 'numpy_solver')
    assert_array_almost_equal(a, sol_cholesky)

    a = CHOLESKY(A, B, 'numpy')
    assert_array_almost_equal(a, sol_cholesky)

    a = CHOLESKY(A, B, 'scipy')
    assert_array_almost_equal(a, sol_cholesky)

    try:
        a = CHOLESKY(A, B, 'dummyy')
        assert False
    except:
        assert True
