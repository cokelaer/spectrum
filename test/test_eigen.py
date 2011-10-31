from spectrum import eigen
import numpy
#from spectrum import MINEIGVAL
from eigen import *
from nose.tools import assert_almost_equal



def test_mineigval():
    tol = 1e-10
    T0=3
    T = numpy.array([-2+.5j, .7-1j],dtype=complex)
    eigval, eigvec = MINEIGVAL(T0 , T, tol)
    print 'Eigenvalue=',eigval
    print 'Eigenvector=',eigvec

    assert_almost_equal(eigval, .488694078106)

    expected_eigvec = numpy.array([ 0.13790622 -1.74155903e-02j , 0.21272177 -4.65701963e-18j,  0.13790622 +1.74155903e-02j])

    assert_almost_equal(eigvec.all(), expected_eigvec.all())


