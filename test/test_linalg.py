


import numpy
from spectrum import linalg
from linalg import * 
from numpy.testing import *

def test_pascal():
    res = pascal(6) 

    assert list(res[5]) == list([   1.,    6.,   21.,   56.,  126.,  252.])
    try:
        pascal(0)
        assert False
    except:
        assert True


def test_svd():
    a = numpy.matrix([[2.8-.4j, 1.6],[3.6-1.2j, 2.4-1.8j],[2+.4j, 2.8-4.4j]])
    U, S, V = csvd(a)
    assert_array_almost_equal(S, numpy.array([ 7.51711296,  2.96867189]))

    csvd(a)


def test_corrmtx():
    C = corrmtx([1,2,3,4,5,6,7,8+1j], 2, method='modified')
    assert_array_almost_equal, C, numpy.array([[ 3.+0.j,  2.+0.j,  1.+0.j],
       [ 4.+0.j,  3.+0.j,  2.+0.j],
       [ 5.+0.j,  4.+0.j,  3.+0.j],
       [ 6.+0.j,  5.+0.j,  4.+0.j],
       [ 7.+0.j,  6.+0.j,  5.+0.j],
       [ 8.+1.j,  7.+0.j,  6.+0.j],
       [ 1.-0.j,  2.-0.j,  3.-0.j],
       [ 2.-0.j,  3.-0.j,  4.-0.j],
       [ 3.-0.j,  4.-0.j,  5.-0.j],
       [ 4.-0.j,  5.-0.j,  6.-0.j],
       [ 5.-0.j,  6.-0.j,  7.-0.j],
       [ 6.-0.j,  7.-0.j,  8.-1.j]])

    C = corrmtx([1,2,3,4,5,6,7,8+1j], 3, method='autocorrelation')
    assert_array_almost_equal(C, numpy.array([[ 1.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 2.+0.j,  1.+0.j,  0.+0.j,  0.+0.j],
       [ 3.+0.j,  2.+0.j,  1.+0.j,  0.+0.j],
       [ 4.+0.j,  3.+0.j,  2.+0.j,  1.+0.j],
       [ 5.+0.j,  4.+0.j,  3.+0.j,  2.+0.j],
       [ 6.+0.j,  5.+0.j,  4.+0.j,  3.+0.j],
       [ 7.+0.j,  6.+0.j,  5.+0.j,  4.+0.j],
       [ 8.+1.j,  7.+0.j,  6.+0.j,  5.+0.j],
       [ 0.+0.j,  8.+1.j,  7.+0.j,  6.+0.j],
       [ 0.+0.j,  0.+0.j,  8.+1.j,  7.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j,  8.+1.j]]))


    C  = corrmtx([1,2,3,4,5,6,7,8+1j], 3, method='covariance')
    assert_array_almost_equal(C, numpy.array([[ 4.+0.j,  3.+0.j,  2.+0.j,  1.+0.j],
       [ 5.+0.j,  4.+0.j,  3.+0.j,  2.+0.j],
       [ 6.+0.j,  5.+0.j,  4.+0.j,  3.+0.j],
       [ 7.+0.j,  6.+0.j,  5.+0.j,  4.+0.j],
       [ 8.+1.j,  7.+0.j,  6.+0.j,  5.+0.j]]))

    C  = corrmtx([1,2,3,4,5,6,7,8+1j], 3, method='postwindowed')
    assert_array_almost_equal(C, numpy.array([[ 4.+0.j,  3.+0.j,  2.+0.j,  1.+0.j],
       [ 5.+0.j,  4.+0.j,  3.+0.j,  2.+0.j],
       [ 6.+0.j,  5.+0.j,  4.+0.j,  3.+0.j],
       [ 7.+0.j,  6.+0.j,  5.+0.j,  4.+0.j],
       [ 8.+1.j,  7.+0.j,  6.+0.j,  5.+0.j],
       [ 0.+0.j,  8.+1.j,  7.+0.j,  6.+0.j],
       [ 0.+0.j,  0.+0.j,  8.+1.j,  7.+0.j],
       [ 0.+0.j,  0.+0.j,  0.+0.j,  8.+1.j]]))

    C  = corrmtx([1,2,3,4,5,6,7,8+1j], 3, method='prewindowed')
    assert_array_almost_equal(C, numpy.array([[ 1.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],
       [ 2.+0.j,  1.+0.j,  0.+0.j,  0.+0.j],
       [ 3.+0.j,  2.+0.j,  1.+0.j,  0.+0.j],
       [ 4.+0.j,  3.+0.j,  2.+0.j,  1.+0.j],
       [ 5.+0.j,  4.+0.j,  3.+0.j,  2.+0.j],
       [ 6.+0.j,  5.+0.j,  4.+0.j,  3.+0.j],
       [ 7.+0.j,  6.+0.j,  5.+0.j,  4.+0.j],
       [ 8.+1.j,  7.+0.j,  6.+0.j,  5.+0.j]]))


    for method in ['prewindowed','postwindowed','autocorrelation','modified','covariance']:
        corrmtx([1,2,3,4,5,6,7,8], 3, method=method)


    try:
        corrmtx([1,2,3,4,5,6,7,8], 3, method='dummy')
        assert False
    except: 
        assert True
